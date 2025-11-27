setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
library(WGCNA)
library(reshape2)
library(tidyverse)
library(stringr)

TISSUE_A <- 'pit'
TISSUE_B <- 'gonad'
SEX <- 'F'
TISSUE_A_RESULTS <- 'DataChapterThree/QENIE/Results/MalesAndFemales/pit_to_gonad.tsv'
TISSUE_B_RESULTS <- 'DataChapterThree/QENIE/Results/Females/gonad_to_pit.tsv'
TISSUE_A_FILE <- "DataChapterThree/Data/salmon.merged_gene_counts_%TISSUE_length_scaled.tsv"
TISSUE_B_FILE <- "DataChapterOne/Data/salmon.merged_gene_counts_%TISSUE_length_scaled.tsv"
SECRETED_PROTEINS <- "DataChapterThree/QENIE/secreted_proteins.txt"
TISSUE_DEGS <- "DataChapterThree/DESeq/DEG_lists/%TISSUE %SEX.txt"
FILTER_FRACTION <- 1/2  # Minimum fraction of non-missing samples for a gene to be considered good.
MIN_EXPRESSION <- 0.9  # At what point a gene is considered not expressed
GENES_IN_GENOME <- 16581


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

n_genes_a <- ncol(tissue_a)
n_genes_b <- ncol(tissue_b)

results_a <- read.csv(TISSUE_A_RESULTS, sep='\t')
results_b <- read.csv(TISSUE_B_RESULTS, sep='\t')

normalise_ssec <- function(ssec, n_genes_tissue, n_genes_genome) {
  return((ssec / n_genes_tissue) * n_genes_genome)
}

results_a$Ssec_norm <- apply(results_a[c(2)], c(1), normalise_ssec, n_genes_b, GENES_IN_GENOME)
results_b$Ssec_norm <- apply(results_b[c(2)], c(1), normalise_ssec, n_genes_a, GENES_IN_GENOME)

write.table(results_a, TISSUE_A_RESULTS, sep = '\t', row.names=FALSE)
write.table(results_b, TISSUE_B_RESULTS, sep = '\t', row.names=FALSE)
