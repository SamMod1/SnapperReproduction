setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
library(WGCNA)
library(reshape2)
library(tidyverse)
library(stringr)

TISSUE_A <- 'pit'
TISSUE_B <- 'gonad'
SEX <- 'all'
TISSUE_A_FILE <- "DataChapterThree/Data/salmon.merged_gene_counts_%TISSUE_length_scaled.tsv"
TISSUE_B_FILE <- "DataChapterOne/Data/salmon.merged_gene_counts_%TISSUE_length_scaled.tsv"
TISSUE_DEGS <- "DataChapterThree/DESeq/DEG_lists/%TISSUE %SEX.txt"
SECRETED_PROTEINS <- "DataChapterThree/QENIE/secreted_proteins.txt"
FILTER_FRACTION <- 1/2  # Minimum fraction of non-missing samples for a gene to be considered good.
MIN_EXPRESSION <- 0.9  # At what point a gene is considered not expressed
ALPHA <- 0.05


tissue_a_file <- str_replace(TISSUE_A_FILE, '%TISSUE', TISSUE_A)
tissue_b_file <- str_replace(TISSUE_B_FILE, '%TISSUE', TISSUE_B)
tissue_a_degs <- str_replace(TISSUE_DEGS, '%TISSUE', TISSUE_A) %>% str_replace('%SEX', SEX)
tissue_b_degs <- str_replace(TISSUE_DEGS, '%TISSUE', TISSUE_B) %>% str_replace('%SEX', SEX)
tissue_a_degs <- scan(tissue_a_degs, what=character())
tissue_b_degs <- scan(tissue_b_degs, what=character())

prep_counts <- function(counts) {
  counts_table <- as.data.frame(t(counts[, -c(1, 2)]))
  names(counts_table) <- counts$gene_name
  rownames(counts_table) <- names(counts[, -c(1, 2)])
  
  if (!SEX == 'all') {
    counts_table <- counts_table[grepl(SEX, row.names(counts_table)),]
  }
  
  return(counts_table)
}


# Import secreted peptides
secreted_proteins <- scan(SECRETED_PROTEINS, what=character())

# Import your data
tissue_a <- read.delim(tissue_a_file, check.names=F)
tissue_b <- read.delim(tissue_b_file, check.names=F)
tissue_a <- prep_counts(tissue_a)
tissue_b <- prep_counts(tissue_b)
tissue_a <- tissue_a[row.names(tissue_a) %in% row.names(tissue_b),]
tissue_b <- tissue_b[row.names(tissue_b) %in% row.names(tissue_a),]

print(head(tissue_a[,colnames(tissue_a)[1:10]]))
print(head(tissue_b[,colnames(tissue_b)[1:10]]))

good_genes_1 = goodSamplesGenes(
  tissue_a, 
  verbose = 3, 
  minFraction = FILTER_FRACTION, 
  tol = MIN_EXPRESSION
)

good_genes_2 = goodSamplesGenes(
  tissue_b, 
  verbose = 3, 
  minFraction = FILTER_FRACTION, 
  tol = MIN_EXPRESSION
)

#counts_table = counts_table[good_genes$goodSamples, good_genes$goodGenes]
bad_genes_overlap <- !(good_genes_1$goodGenes) & !(good_genes_2$goodGenes)

tissue_a <- tissue_a[!bad_genes_overlap]
tissue_b <- tissue_b[!bad_genes_overlap]

if (!all(good_genes_1$goodSamples) | !all(good_genes_1$goodSamples)) {
  print('The following sample(s) were flagged as unsuitable by WGCNA:')
  print(rownames(counts_table)[!good_genes$goodSamples])
  stop('Sample(s) were flagged as unsuitable by WGCNA')
}

secreted_proteins <- secreted_proteins[secreted_proteins %in% colnames(tissue_a)]

# genes <- unique(append(sample(colnames(tissue_a), 1000), secreted_proteins))#sample(secreted_proteins, 500)))
# tissue_a <- tissue_a[genes]
# tissue_b <- tissue_b[genes[genes %in% colnames(tissue_b)]]

# Using these two datasets we will proceed:

# Construct cross-tissue correlation and pvalue matrices
correlations = bicorAndPvalue(tissue_a, tissue_b, use='pairwise.complete.obs')

# Compute rowsum -log(pvalue) for ranking interactions and filter for factors proteins as secreted
scores_tissue_1 = rowSums(-log(correlations$p), na.rm = TRUE)
scores_tissue_2 = colSums(-log(correlations$p), na.rm = TRUE)

# Retain only secreted peptides, normalize the Ssec score by number of target tissue probes (In this case, adipose contains 12242 genes)
# order by "sig score"
score_df_1 = data.frame(Gene_symbol = names(scores_tissue_1), score = scores_tissue_1) %>%
  filter(Gene_symbol %in% secreted_proteins) %>%
  mutate(score_df_1 = score / length(colnames(tissue_a))) %>%
  select(-score) %>%
  arrange(desc(score_df_1))

score_df_2 = data.frame(Gene_symbol = names(scores_tissue_2), score = scores_tissue_2) %>%
  filter(Gene_symbol %in% secreted_proteins) %>%
  mutate(score_df_2 = score / length(colnames(tissue_a))) %>%
  select(-score) %>%
  arrange(desc(score_df_2))

#This produces a table to each secreted protein and its respective significance score across adipose transcripts
#Note that Notum is listed as the 5th with an Sssec of 4.148
#For our pipeline, we next check the tissue-specificty using BioGPS - note that this step is not necessary, but makes us more confident when conditioning the pathway enrichment.  This moves Notum up to the 2nd ranked protein

#Condition correlation matrix on a by-gene basis for pathway enrichment - this example will focus on the protein, Notum
#remove gene of interest (Notum) from the correlation matrix
#the rownames of the file correspond to gene symbols for pathway enrichment, whereas the second column contains the bicor coefficent

#saveRDS(correlations, file='corr.RDS')
bicor = data.frame(melt(correlations$bicor))
p <- data.frame(melt(correlations$p))
rm(correlations)
colnames(bicor) = c('Gene_1', 'Gene_2', 'bicor')
bicor$p <- p$value
#bicor$p_adj <- p.adjust(bicor$p, method = 'BH')

#name <- paste(TISSUE_A, TISSUE_B, sep= '_vs_') %>% paste('.csv', sep='')
bicor <- bicor[(bicor$Gene_1 %in% secreted_proteins) | (bicor$Gene_2 %in% secreted_proteins),]

n_cor <- integer()
for (protein in score_df_1$Gene_symbol) {
  df <- filter(bicor, Gene_1 == protein)
  p_adj <- p.adjust(df$p, method = 'BH')
  n_cor <- append(n_cor, sum(p_adj <= ALPHA, na.rm=TRUE))
}
score_df_1$n_cor <- n_cor

n_cor <- integer()
for (protein in score_df_2$Gene_symbol) {
  df <- filter(bicor, Gene_2 == protein)
  p_adj <- p.adjust(df$p, method = 'BH')
  n_cor <- append(n_cor, sum(p_adj <= ALPHA, na.rm=TRUE))
}
score_df_2$n_cor <- n_cor

score_df_1$is_deg <- score_df_1$Gene_symbol %in% tissue_a_degs
score_df_2$is_deg <- score_df_2$Gene_symbol %in% tissue_b_degs

write.csv(bicor, name, row.names=FALSE)

name_1 <- paste(TISSUE_A, TISSUE_B, sep= '_to_') %>% paste('.tsv', sep='')
name_2 <- paste(TISSUE_B, TISSUE_A, sep= '_to_') %>% paste('.tsv', sep='')
write.table(score_df_1, file=name_1, row.names=F, col.names=T, sep='\t', quote=F)
write.table(score_df_2, file=name_2, row.names=F, col.names=T, sep='\t', quote=F)
