if (!require("pacman")) {install.packages("pacman")}
pacman::p_load(BiocManager, gridExtra, ggplot2, ggdendro, stringr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
if (!('WGCNA' %in% installed.packages())){BiocManager::install("WGCNA")}
library(DESeq2)
library(WGCNA)
options(stringsAsFactors = FALSE)
set.seed(1)
source('functions.R') # Includes:
# binarise_phenotype_data(phenotype_data): Binarises phenotype factors into a format WGCNA can accept
# plot_membership_trait_significance(module, trait, module_membership, trait_significance)
# CorLevelPlot: Modified source code for the CorLevelPlot function, accessed on 10/01/2022 from:
# https://rdrr.io/github/kevinblighe/CorLevelPlot/src/R/CorLevelPlot

#source("C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/GoEnrichment/go_enrichment.R") # Takes a little time to initialise but adds the function:
# enrichment(genes, database)
# enrich_from_module(module, module_membership, database)

#enableWGCNAThreads()

OUTPUT <- 'DataChapterFour/WGCNA/Results'
NAME <- 'gene_counts_gonad'
TISSUE <- 'gonad'
NETWORK_TYPE <- 'signed'
#DIFFERENTIAL_EXPRESSION_DATA_DIR <- "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq2/Output - Sex V Month Wald/DEGs"
COUNTS_FILE <- 'DataChapterTwo/Data/salmon.merged.gene_counts_gonad.tsv'
FILTER_FRACTION <- 1/2  # Minimum fraction of non-missing samples for a gene to be considered good.
MIN_EXPRESSION <- 0.9  # At what point a gene is considered not expressed
BLOCKSIZE <- 27000  # Number of genes to do at once. If below total number of genes results may be effected. 16Gb can do 20,000
ALPHA <- 0.05
LIMIT_SEX <- 'Juvenile'


# =========================== Prep data ===========================

counts <- read.table(COUNTS_FILE, sep='\t', header=TRUE)

counts_table <- as.data.frame(t(counts[, -c(1, 2)]))
names(counts_table) <- counts$gene_name
rownames(counts_table) <- names(counts[, -c(1, 2)])

good_genes = goodSamplesGenes(
  counts_table, 
  verbose = 3, 
  minFraction = FILTER_FRACTION, 
  tol = MIN_EXPRESSION
)

if (!good_genes$allOK) {counts_table = counts_table[good_genes$goodSamples, good_genes$goodGenes]}
if (!all(good_genes$goodSamples)) {
  print('The following sample(s) were flagged as unsuitable by WGCNA:')
  print(rownames(counts_table)[!good_genes$goodSamples])
  stop('Sample(s) were flagged as unsuitable by WGCNA')
}

traits <- pca_dendro_stage(counts_table)

SAMPLES_TO_EXCLUDE <- c('X13_JM')

if (is.character(LIMIT_SEX)) {
  samples_to_exclude <- append(
    SAMPLES_TO_EXCLUDE, 
    row.names(subset(traits, traits[,LIMIT_SEX] == 0))
  )
} else {samples_to_exclude <- SAMPLES_TO_EXCLUDE}

counts_table <- counts_table[!rownames(counts_table) %in% samples_to_exclude,]

traits <- pca_dendro_stage(counts_table)

if (is.character(LIMIT_SEX)) {traits[LIMIT_SEX] <- NULL}

# =========================== Run WGCNA ===========================

#counts_table <- counts_table[,sample(colnames(counts_table), 10000)]
counts_table <- log2(counts_table + 1)
#counts_table <- log(counts_table + 1)
#counts_table <- varianceStabilizingTransformation(counts_table)

powers = c(c(1:10), seq(from = 12, to=30, by=2))

soft_thresholds = pickSoftThreshold(
  counts_table, 
  powerVector = powers, 
  verbose = 5, 
  networkType = NETWORK_TYPE,
  blockSize=BLOCKSIZE / 3
)

ggplot(
  data.frame('x' = soft_thresholds$fitIndices[,1], 'y' = -sign(soft_thresholds$fitIndices[,3])*soft_thresholds$fitIndices[,2]),
  aes(x = x, y = y, label = x)
) + 
  geom_text(color = 'red') + ylab('Scale free topology model fit, signed R^2') + xlab('Power') +
  geom_hline(yintercept = 0.90, col = 'green') + geom_hline(yintercept = 0.80, col = 'orange') +
  geom_hline(yintercept = 0.70, col = 'red')

ggplot(
  data.frame('x' = soft_thresholds$fitIndices[,1], 'y' = soft_thresholds$fitIndices[,5]),
  aes(x = x, y = y, label = x)
) + geom_text(color = 'red') + ylab('Mean Connectivity') + xlab('Power')

while (TRUE) {
  #soft_threshold <- readline(prompt="Choose soft threshold: ")
  tryCatch(
    {
      soft_threshold <- 16 #as.numeric(soft_threshold)
      break
    }
    , warning = function(e) {print('Soft threshold must be numeric!')}
  )
}

network <- blockwiseModules(counts_table, power = soft_threshold,
                       TOMType = NETWORK_TYPE, minModuleSize = 50,
                       mergeCutHeight = 0.25, #detectCutHeight = 0.99,
                       maxPOutliers = 0.05, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, randomSeed = 1,
                       maxBlockSize = BLOCKSIZE,
                       saveTOMFileBase = NAME, 
                       verbose = 3, corType = "bicor")

merged_colors = network$colors
plotDendroAndColors(network$dendrograms[[1]], merged_colors[network$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

module_colors = network$colors
eigengenes = network$MEs
gene_tree = network$dendrograms[[1]];
save(eigengenes, module_colors, gene_tree, counts_table, traits,
     file = paste(NAME, ".RData", sep=''))

saveRDS(network, paste(OUTPUT, 'wgcna_network.rds', sep='/'))

write.csv(data.frame(module_colors), paste(OUTPUT, 'module_membership.csv', sep='/'))
write.csv(eigengenes, paste(OUTPUT, 'eigengenes.csv', sep='/'))


# =========================== Correlation with traits ===========================

#network <- readRDS(paste(OUTPUT, 'wgcna_network.rds', sep='/'))
#names <- load(file=paste(OUTPUT, paste(NAME, '.RData', sep=''), sep='/'))

n_genes = ncol(counts_table);
n_samples = nrow(counts_table);
# Recalculate MEs with color labels
eigengenes0 = moduleEigengenes(counts_table, module_colors)$eigengenes
eigengenes = orderMEs(eigengenes0)
trait_correlations = cor(eigengenes, traits, use = "p")
trait_pvalue = corPvalueStudent(trait_correlations, n_samples)

heatmap_data <- merge(eigengenes, traits, by = 'row.names')

CorLevelPlot(heatmap_data,
             x = colnames(trait_pvalue),
             y = colnames(eigengenes),
             col = c('steelblue', 'skyblue', 'white', 'pink', 'brown2'),
             yLabelColors = substring(colnames(eigengenes), 3),
             titleX = 'Condition')

module_membership_cor = as.data.frame(cor(counts_table, eigengenes, use = "p"))
module_membership_p = as.data.frame(corPvalueStudent(as.matrix(module_membership_cor), n_samples))
colnames(module_membership_p) <- paste(colnames(module_membership_p), 'p', sep='_')
module_membership <- cbind(module_membership_cor, module_membership_p)
col_order <- c()
for (module in colnames(module_membership_cor)) {
  col_order <- append(col_order, module)
  col_order <- append(col_order, paste(module, 'p', sep='_'))
}
module_membership <- module_membership[,col_order]

trait_significance <- data.frame(row.names = colnames(counts_table))
for (trait in colnames(traits)) {
  trait_data <- traits[trait]
  
  correlation <- cor(counts_table, trait_data, use = "p")
  significance <- corPvalueStudent(as.matrix(correlation), n_samples)
  colnames(significance) <- paste(trait, 'p', sep='_')
  trait_significance <- cbind(trait_significance, correlation, significance)
}

module_colors_df <- read.csv(paste(OUTPUT, 'module_membership.csv', sep='/'), row.names=1)
module_membership <- cbind(module_membership, module_colors_df)

head(module_membership)
head(trait_significance)

write.csv(trait_pvalue, paste(OUTPUT, 'eigengene_trait_significance.csv', sep='/'))
write.csv(module_membership, paste(OUTPUT, 'module_membership_vals.csv', sep='/'))
write.csv(trait_significance, paste(OUTPUT, 'trait_significance.csv', sep='/'))
hubs <- chooseTopHubInEachModule(counts_table, module_colors)

save.image(file=paste(OUTPUT, 'whole_analysis.RData', sep='/'))


# =========================== Filter module genes by differential expression ===========================


load(paste(OUTPUT, 'whole_analysis.RData', sep='/'))


module <- 'blue'

go <- enrich_from_module(module, module_membership, 'GO', dotplot=FALSE, alpha=ALPHA, gene_universe=colnames(counts_table))
kegg1 <- enrich_from_module(module, module_membership, 'KEGG', dotplot=FALSE, alpha=ALPHA, gene_universe=colnames(counts_table))
kegg <- kegg_from_module(module, module_membership, pathway=TRUE, dotplot=FALSE, alpha=ALPHA)


# =========================== Final plotting and analysis ===========================

modules <- c()
n_genes <- c()
top_go <- c()
top_kegg <- c()
go_readable <- c()
kegg_readable <- c()
idx <- 1

for (module in unique(module_membership$module_colors)) {
  modules[idx] <- module
  genes <- rownames(subset(module_membership, module_membership$module_colors == module))
  go <- enrichment(genes, 'GO', dotplot=FALSE, alpha=ALPHA, gene_universe=colnames(counts_table))
  kegg <- kegg_enrichment(genes, pathway=TRUE, dotplot=FALSE, alpha=ALPHA, gene_universe=colnames(counts_table))
  
  n_genes[idx] <- length(genes)
  
  if (!is.null(go)) {
    go <- subset(go@result, go@result$p.adjust <= 0.05)
    go <- go[1:10,]
    go <- go[!is.na(go$ID),]
    go_readable[idx] <- paste(as.character(go$Description), collapse = '; ')
    top_go[idx] <- paste(as.character(go$ID), collapse = '; ')
  } else {top_go[idx] <- NA; go_readable[idx] <- NA}
  
  if (!is.null(kegg)) {
    kegg <- subset(kegg@result, kegg@result$p.adjust <= 0.05)
    kegg <- kegg[1:10,]
    kegg <- kegg[!is.na(kegg$ID),]
    kegg_readable[idx] <- paste(as.character(kegg$Description), collapse = '; ')
    top_kegg[idx] <- paste(as.character(kegg$ID), collapse = '; ')
  } else {top_kegg[idx] <- NA; kegg_readable[idx] <- NA}
  
  idx <- idx + 1
}

df <- data.frame(
  module = modules,
  n_genes = n_genes,
  top_go = top_go,
  top_kegg = top_kegg,
  go_readable = go_readable,
  kegg_readable = kegg_readable
)

write.csv(df, paste(OUTPUT, 'functional_analysis.csv', sep='/'), row.names = FALSE)

# =========================== Reload previous analysis ===========================

if (!require("pacman")) {install.packages("pacman")}
pacman::p_load(BiocManager, gridExtra, ggplot2, ggdendro, stringr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
if (!('WGCNA' %in% installed.packages())){BiocManager::install("WGCNA")}
library(WGCNA)
options(stringsAsFactors = FALSE)
set.seed(1)
source('functions.R') # Includes:
# binarise_phenotype_data(phenotype_data): Binarises phenotype factors into a format WGCNA can accept
# plot_membership_trait_significance(module, trait, module_membership, trait_significance)
# CorLevelPlot: Modified source code for the CorLevelPlot function, accessed on 10/01/2022 from:
# https://rdrr.io/github/kevinblighe/CorLevelPlot/src/R/CorLevelPlot

#source("C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/GoEnrichment/go_enrichment.R") # Takes a little time to initialise but adds the function:
# enrichment(genes, database)
# enrich_from_module(module, module_membership, database)

#enableWGCNAThreads()

OUTPUT <- 'DataChapterFour/WGCNA/Results'
NAME <- 'gene_counts_gonad'
TISSUE <- 'gonad'
NETWORK_TYPE <- 'signed'
#DIFFERENTIAL_EXPRESSION_DATA_DIR <- "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq2/Output - Sex V Month Wald/DEGs"
COUNTS_FILE <- 'DataChapterTwo/Data/salmon.merged.gene_counts_gonad.tsv'
FILTER_FRACTION <- 1/2  # Minimum fraction of non-missing samples for a gene to be considered good.
MIN_EXPRESSION <- 0.9  # At what point a gene is considered not expressed
BLOCKSIZE <- 27000  # Number of genes to do at once. If below total number of genes results may be effected. 16Gb can do 20,000
ALPHA <- 0.05
LIMIT_SEX <- 'Juvenile'

# =========================== Prep data ===========================

counts <- read.table(COUNTS_FILE, sep='\t', header=TRUE)

counts_table <- as.data.frame(t(counts[, -c(1, 2)]))
names(counts_table) <- counts$gene_name
rownames(counts_table) <- names(counts[, -c(1, 2)])

good_genes = goodSamplesGenes(
  counts_table, 
  verbose = 3, 
  minFraction = FILTER_FRACTION, 
  tol = MIN_EXPRESSION
)

if (!good_genes$allOK) {counts_table = counts_table[good_genes$goodSamples, good_genes$goodGenes]}
if (!all(good_genes$goodSamples)) {
  print('The following sample(s) were flagged as unsuitable by WGCNA:')
  print(rownames(counts_table)[!good_genes$goodSamples])
  stop('Sample(s) were flagged as unsuitable by WGCNA')
}

traits <- pca_dendro_stage(counts_table)

SAMPLES_TO_EXCLUDE <- c('X13_JM')

if (is.character(LIMIT_SEX)) {
  samples_to_exclude <- append(
    SAMPLES_TO_EXCLUDE, 
    row.names(subset(traits, traits[,LIMIT_SEX] == 0))
  )
} else {samples_to_exclude <- SAMPLES_TO_EXCLUDE}

counts_table <- counts_table[!rownames(counts_table) %in% samples_to_exclude,]

traits <- pca_dendro_stage(counts_table)

if (is.character(LIMIT_SEX)) {traits[LIMIT_SEX] <- NULL}

load(paste(OUTPUT, 'whole_analysis.RData', sep='/'))

n_genes = ncol(counts_table);
n_samples = nrow(counts_table);
# Recalculate MEs with color labels
eigengenes0 = moduleEigengenes(counts_table, module_colors)$eigengenes
eigengenes = orderMEs(eigengenes0)
trait_correlations = cor(eigengenes, traits, use = "p")
trait_pvalue = corPvalueStudent(trait_correlations, n_samples)

heatmap_data <- merge(eigengenes, traits, by = 'row.names')

CorLevelPlot(heatmap_data,
             x = colnames(trait_pvalue), #c('F1', 'F2', 'F3', 'F4', 'F5'),#
             y = colnames(eigengenes),
             col = c('steelblue', 'skyblue', 'white', 'pink', 'brown2'),
             yLabelColors = substring(colnames(eigengenes), 3),
             titleX = 'Condition')

hubs <- chooseTopHubInEachModule(counts_table, module_colors)
# trait_name <- 'F4'
# 
# plot_membership_trait_significance(module, trait_name, module_membership, trait_significance)
