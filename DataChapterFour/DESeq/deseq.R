setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
if (!require("pacman")) install.packages("pacman")
pacman::p_load(BiocManager, ggplot2, edgeR, stringr)
if (!('DESeq2' %in% installed.packages())){BiocManager::install("DESeq2")}
library(DESeq2)
source('functions.R')
set.seed(1)

ALPHA <- 0.05
EXPRESSION_CUTOFF <- 1
EXPRESSION_CUTOFF_FRACTION <- 1/2 # Fraction of samples with no expression to filter out
COUNTS_FILE <- 'DataChapterTwo/Data/salmon.merged.gene_counts_gonad.tsv'
OUTPUT_DIR <- 'DataChapterFour/DESeq/Results'
SAMPLE_METADATA <- 'DataChapterOne/Data/sampling_data_juveniles.csv'
SAMPLES_TO_EXCLUDE = c('X90_M3', 'X13_M1')
PLOT_WIDTH = 1200
PLOT_HEIGHT = 800

# ===================Init data=============================================

origional_counts <- read.table(COUNTS_FILE, sep='\t', header=TRUE)
sample_metadata <- read.csv(SAMPLE_METADATA)

samples <- colnames(origional_counts[, 3:ncol(origional_counts)])
sample_metadata <- extract_sample_metadata(samples, sample_metadata)
sample_metadata <- subset(sample_metadata, !sample_metadata$sample %in% SAMPLES_TO_EXCLUDE)
sample_metadata <- subset(sample_metadata, sample_metadata$sex == 'JM' | sample_metadata$sex == 'JF')
adults <- sample_metadata$sample
origional_counts <- origional_counts[colnames(origional_counts) %in% append(c('gene_id', 'gene_name'), adults)]

row.names(origional_counts) <- origional_counts$gene_id
origional_counts <- origional_counts[!colnames(origional_counts) %in% c('gene_id', 'gene_name')]
for (col in colnames(origional_counts)){origional_counts[,col] <- as.integer(round(origional_counts[,col]))}
row.names(sample_metadata) <- sample_metadata$sample
sample_metadata <- sample_metadata[,2:ncol(sample_metadata)]
sample_metadata$sex <- as.factor(sample_metadata$sex)
sample_metadata$stage <- as.factor(sample_metadata$stage)

# ===========================Preprocessing===================================

counts <- DGEList(counts=origional_counts)

genes_to_keep <- filterByExpr(counts)

total_genes <- nrow(origional_counts)
print('Genes filtered:')
print(total_genes - sum(genes_to_keep))
print('Number of genes after filering:')
print(sum(genes_to_keep))

normalisation_dataset <- DESeqDataSetFromMatrix(
  countData=origional_counts,
  colData=sample_metadata,
  design=~stage
)
normalisation_dataset <- estimateSizeFactors(normalisation_dataset)
normalised_counts <- counts(normalisation_dataset, normalized=TRUE)
write.table(normalised_counts, 'normalised_counts.tsv', sep = '\t', row.names=FALSE)

# ======================sex analysis=========================================

sample_sex <- factor(t(sample_metadata["sex"][[1]]))

sex_counts <- origional_counts

sex_dataset <- DESeqDataSetFromMatrix(
  countData=sex_counts,
  colData=sample_metadata,
  design=~sex
)

sex_dataset <- estimateSizeFactors(sex_dataset)
sex_dataset <- sex_dataset[genes_to_keep,]

pca_plt <- create_pca_plot(sex_dataset, "sex", "gonad")
ggsave(paste(OUTPUT_DIR, 'Sex PCA.png', sep = '/'), plot=pca_plt, width = PLOT_WIDTH * 2, height = PLOT_HEIGHT * 2, units='px')

wald_dataset <- DESeq(sex_dataset, test='Wald')

outcome <- results(wald_dataset, test="Wald")

print(summary(outcome))
test_results <- subset(outcome, outcome$padj < ALPHA)
test_results <- test_results[order(test_results$padj),]  
print(head(test_results))

test_results <- as.data.frame(test_results)
test_results$gene <- row.names(test_results)
print('Number of DEGs:')
print(sum(test_results$padj <= ALPHA))
write.table(test_results, paste(OUTPUT_DIR, 'Sex Wald Results.tsv', sep = '/'), sep='\t', row.names=FALSE)
