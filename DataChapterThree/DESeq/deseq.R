if (!require("pacman")) install.packages("pacman")
pacman::p_load(BiocManager, ggplot2, edgeR, stringr)
if (!('DESeq2' %in% installed.packages())){BiocManager::install("DESeq2")}
library(DESeq2)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
set.seed(1)


ALPHA <- 0.05
LOG_CHANGE_CUTOFF <- 1
EXPRESSION_CUTOFF <- 1 # If the expression cutoff fraction number of samples is below this the gene is removed
EXPRESSION_CUTOFF_FRACTION <- 1/2 # Fraction of samples with no expression to filter out
COUNTS_FILES <- '/DataChapterThree/Data/salmon.merged.gene_counts_%tissue.tsv'
OUTPUT_DIR <- 'ResultsStage'
SAMPLE_METADATA <- 'DataChapterOne/Data/sampling_data.csv'

ANALYSES <- list(list('gonad'), list('pit', 'mid', 'brain'))
FACTORS <- list(list('sex', 'stage', 'stage'),
           list('sex', 'stage', 'stage'))
FULL_DESIGN <- sex~stage
DESIGNS <- list(list(~sex, ~stage, ~stage),
                           list(~sex, ~stage+month, ~stage+month))
INCLUDE_FACTOR_LEVELS <- list(list(list('F', 'M', 'J'), list('F1', 'F2', 'F3', 'F4', 'F5'), list('M2', 'M3', 'M4')),
                        list(list('F', 'M'), list('F1', 'F3', 'F4', 'F5'), list('M3', 'M4')))


SAMPLES_TO_EXCLUDE = c('X90_M3', 'X13_M1')

PLOT_WIDTH = 1200
PLOT_HEIGHT = 800


volcano_plot <- function(log2_fold, p_vals, alpha, log_change_cutoff, title) {
  
  df <- data.frame(log2_fold <- log2_fold, p_val <- p_vals)
  df$result <- 'No Change'
  df$result[df$p_val < alpha & df$log2_fold > log_change_cutoff] <- 'Upregulated'
  df$result[df$p_val < alpha & df$log2_fold < -log_change_cutoff] <- 'Downregulated'
  fig <- ggplot(df, aes(x=log2_fold, y=-log10(p_val), color=result)) + 
    geom_point() + scale_color_manual(values=c('blue', 'lightblue4', 'red')) +
    ggtitle(title)
  
  return(fig)
}


create_pca_plot <- function(dds, factor, tissue) {
  vsd <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup=c(factor), returnData=TRUE)
  colnames(pca_data)[colnames(pca_data) == factor] <- 'variable'
  pca_data$variable <- as.character(pca_data$variable)
  plt <- ggplot(pca_data, aes(PC1, PC2, color=variable, text=name)) + 
      geom_point() + scale_color_discrete(name = factor) + ggtitle(tissue)
  
  write.csv(pca_data, file='temp.csv')
  
  return(plt)
}


calculate_ps <- function(deseq_dataset, variables, factor, alpha) {
  outcome <- results(deseq_dataset, contrast=c(factor, variables[1], variables[2]))
  
  print(summary(outcome))
  significant_results <- subset(outcome, outcome$padj < alpha)
  significant_results <- significant_results[order(significant_results$padj),]

  if (nrow(significant_results) > 0) {
    significant_results$name <- toString(variables)
  } else {
    significant_results$name <- c()
  }
  
  print(significant_results)
  
  return(list(outcome, significant_results))
}


extract_sample_metadata <- function(samples, sample_metadata) {
  sample_numbers <- as.numeric(sub('X', '', str_extract(samples, regex('X[0-9]*'))))
  metadata <- subset(sample_metadata, sample_metadata$Fish.ID %in% sample_numbers)
  metadata <- metadata[match(sample_numbers, metadata$Fish.ID),]
  
  df <- data.frame(
    sample = samples, 
    sex = metadata$Sex, 
    stage = metadata$Stage, 
    month = str_replace_all(paste(metadata$Date, metadata$Sex), '/', '-')
  )
  
  return(df)
}


preprocessing <- function(counts_df, sample_factors) {
  
  counts <- DGEList(counts=counts_df, group=sample_factors)
  
  genes_to_keep <- filterByExpr(counts)
  
  total_genes <- nrow(counts_df)
  counts_df <- counts_df[genes_to_keep,]
  print('Genes filtered:')
  print(total_genes - nrow(counts_df))
  print('Number of genes after filering:')
  print(nrow(counts_df))
  
  return(counts_df)
}


analyse <- function(
    wald_dataset, 
    lrt_dataset,
    factor, 
    factor_levels, 
    alpha, 
    tissue, 
    output, 
    pca_plt, 
    fold_change_cutoff, 
    expression_cutoff_fraction, 
    expression_cutoff
  ) {
  
  data_dir <- paste(output, 'Data', sep='//')
  if (!file.exists(data_dir)) {
    dir.create(data_dir)
  }
  alpha_str <- paste('Alpha=', as.character(alpha), sep='')
  output <- paste(output, alpha_str, sep='//')
  if (!file.exists(output)) {
    dir.create(output)
  }
  plots_dir <- paste(output, 'Plots', sep='//')
  if (!file.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  volcano_dir <- paste(output, paste('Log Change=', as.character(fold_change_cutoff), sep=''), sep='//')
  if (!file.exists(volcano_dir)) {
    dir.create(volcano_dir)
  }
  
  pca_filename <- paste(paste(plots_dir, tissue, sep='//'), factor)
  pca_filename <- paste(paste(pca_filename, toString(factor_levels)), 'PCA.png')
  ggsave(pca_filename, plot=pca_plt, width = PLOT_WIDTH * 2, height = PLOT_HEIGHT * 2, units='px')
  data_filename <- paste(paste(data_dir, tissue, sep='//'), factor)
  data_filename <- paste(paste(data_filename, toString(factor_levels)), 'tsv', sep='.')
  
  outcome <- results(lrt_dataset, test="LRT")
  
  print(summary(outcome))
  test_results <- subset(outcome, outcome$padj < alpha)
  test_results <- test_results[order(test_results$padj),]  
  print(head(test_results))
  
  test_results <- as.data.frame(test_results)
  
  test_results$gene <- row.names(test_results)
  print('Number of DEGs:')
  print(sum(
    test_results$padj <= alpha &
      test_results$log2FoldChange <= -fold_change_cutoff |
      test_results$log2FoldChange >= fold_change_cutoff
  ))
  write.table(test_results, data_filename, sep='\t', row.names=FALSE)
  
  for (combination in combn(unlist(factor_levels), 2, simplify=FALSE)) {

    ma_filename <- paste(paste(plots_dir, tissue, sep='//'), factor)
    ma_filename <- paste(paste(paste(paste(ma_filename, combination[1]), 'vs'), combination[2]), 'MA.png')
    csv_filename <- paste(paste(data_dir, tissue, sep='//'), factor) 
    csv_filename <- paste(paste(paste(paste(csv_filename, combination[1]), 'vs'), combination[2]), '.tsv', sep='')
    volcano_filename <- paste(paste(volcano_dir, tissue, sep='//'), factor)
    volcano_filename <- paste(paste(paste(paste(volcano_filename, combination[1]), 'vs'), combination[2]), 'volcano.png')
    
    # Filter out low expression:
    factor1 <- unlist(lrt_dataset@colData[factor]) == combination[1]
    factor2 <- unlist(lrt_dataset@colData[factor]) == combination[2]
    cutoff_f1 <- sum(factor1) * expression_cutoff_fraction
    cutoff_f2 <- sum(factor2) * expression_cutoff_fraction
    unfiltered <- nrow(lrt_dataset)
    #cutoff_number <- ncol(lrt_dataset) * expression_cutoff_fraction
    counts_df <- counts(lrt_dataset)
    keep <- c()
    
    for (gene in rownames(counts_df)) {
      row <- counts_df[gene,]
      factor_1_values <- row[factor1]
      factor_2_values <- row[factor2]
      factor_1_above <- factor_1_values >= expression_cutoff
      factor_2_above <- factor_2_values >= expression_cutoff
      factor_1_pass <- sum(factor_1_above) >= cutoff_f1
      factor_2_pass <- sum(factor_2_above) >= cutoff_f2
      keep <- append(keep, c(factor_1_pass | factor_2_pass))
    }
    contrast_dataset <- lrt_dataset[keep,]
    
    filtered <- nrow(contrast_dataset)
    print('Further genes filtered due to low contrast expression:')
    print(unfiltered - filtered)

    p_vals_0 <- calculate_ps(contrast_dataset, as.vector(combination), factor, alpha)
    log_fold <- as.data.frame(p_vals_0[[1]])

    genes_to_keep <- row.names(subset(
      log_fold, 
      log_fold$log2FoldChange <= -fold_change_cutoff |
        log_fold$log2FoldChange >= fold_change_cutoff
    ))
    new_dataset <- contrast_dataset[row.names(contrast_dataset) %in% genes_to_keep,]
    log_fold <- subset(log_fold, !row.names(log_fold) %in% genes_to_keep)
    
    p_vals <- calculate_ps(new_dataset, as.vector(combination), factor, alpha)
    
    results <- as.data.frame(p_vals[[1]])
    results <- base::rbind(results, log_fold)
    results <- results[order(results$padj),]
    print('Number of DEGs:')
    print(sum(
      results$padj <= alpha &
        results$log2FoldChange <= -fold_change_cutoff |
        results$log2FoldChange >= fold_change_cutoff
    ))
    print(head(results))
    
    num_signif_genes <- nrow(subset(results, results$padj < alpha))
    ma_title <- paste(tissue, paste(combination[1], combination[2], sep = ' vs '), sep='\n')
    ma_title <- paste(paste(ma_title, '\nSignificant Genes ='), as.character(num_signif_genes))
    png(ma_filename, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    plotMA(p_vals_0[[1]], main = ma_title)
    dev.off()
    
    results <- results[order(results$padj),]
    
    volcano_title <- paste(combination[1], combination[2], sep=' vs ')
    n_degs <- nrow(subset(results, results$padj < alpha & (
      results$log2FoldChange > fold_change_cutoff | results$log2FoldChange < -fold_change_cutoff)))
    volcano_title <- paste(volcano_title, paste('n-degs =', as.character(n_degs)), sep='\n')
    
    ggsave(volcano_filename, volcano_plot(results$log2FoldChange, results$padj, alpha, fold_change_cutoff, volcano_title))
    
    results$gene <- row.names(results)
    write.table(results, csv_filename, sep='\t', row.names=FALSE)
  }
}


run_analysis <- function(
    counts_df, 
    tissue, 
    factor, 
    design, 
    include_factor_levels, 
    sample_metadata,
    samples_to_exclude, 
    output, 
    alpha, 
    fold_change_cutoff,
    expression_cutoff_fraction,
    expression_cutoff
  ) {
  
  samples <- colnames(counts_df[, 3:ncol(counts_df)])
  sample_metadata <- extract_sample_metadata(samples, sample_metadata)

  if (length(include_factor_levels) > 0) {
   sample_metadata <- subset(sample_metadata, sample_metadata[factor][[1]] %in% include_factor_levels)
  }
  
  sample_metadata <- subset(sample_metadata, !sample_metadata$sample %in% samples_to_exclude)
  
  samples_to_exclude <- append(samples[!samples %in% sample_metadata$sample], samples_to_exclude)
  
  counts_df <- counts_df[!colnames(counts_df) %in% samples_to_exclude]
  
  # Prep counts df for deseq
  row.names(counts_df) <- counts_df$gene_id
  counts_df <- counts_df[!colnames(counts_df) %in% c('gene_id', 'gene_name')]
  for (col in colnames(counts_df)){counts_df[,col] <- as.integer(round(counts_df[,col]))}
  
  row.names(sample_metadata) <- sample_metadata$sample
  sample_metadata <- sample_metadata[,2:ncol(sample_metadata)]
  conditions <- factor(t(sample_metadata[factor][[1]]))
  
  counts_df <- preprocessing(counts_df, conditions)
  
  deseq_dataset <- DESeqDataSetFromMatrix(
    countData=counts_df,
    colData=sample_metadata,
    design=design
  )
  
  deseq_dataset <- estimateSizeFactors(deseq_dataset)
  normalized_counts <- counts(deseq_dataset, normalized=TRUE)
  
  write.table(normalized_counts, 'temp.csv', sep = '\t', row.names=TRUE)
  
  pca_plt <- create_pca_plot(deseq_dataset, factor, tissue)
  
  lrt_dataset <- DESeq(deseq_dataset, test='LRT', reduced = ~1)
  wald_dataset <- DESeq(deseq_dataset, test='Wald')
  
  analyse(
    wald_dataset,
    lrt_dataset,
    factor, 
    include_factor_levels, 
    alpha,
    tissue, 
    output, 
    pca_plt, 
    fold_change_cutoff,
    expression_cutoff,
    expression_cutoff_fraction
  )
}


sample_metadata_full <- read.csv(SAMPLE_METADATA)
for (analysis_idx in 1:length(ANALYSES)) {
  for (tissue_idx in 1:length(ANALYSES[[analysis_idx]])) {
    for (factor_idx in 1:length(FACTORS[[analysis_idx]])) {
      
      
      
      tissue <- ANALYSES[[analysis_idx]][[tissue_idx]]
      factor <- FACTORS[[analysis_idx]][[factor_idx]]
      design <- DESIGNS[[analysis_idx]][[factor_idx]]
      include_factor_levels <- INCLUDE_FACTOR_LEVELS[[analysis_idx]][[factor_idx]]
      
      print('==================')
      print(tissue)
      print(factor)
      print(toString(include_factor_levels))
      print('==================')
      
      counts_file <- gsub('\\%tissue', tissue, COUNTS_FILES)
      counts_df <- read.table(counts_file, sep='\t', header=TRUE)
      
      # full_lrt(
      #   counts_df,
      #   FULL_DESIGN,
      #   OUTPUT_DIR,
      #   ALPHA,
      #   LOG_CHANGE_CUTOFF,
      #   SAMPLES_TO_EXCLUDE
      # )
      run_analysis(
        counts_df,
        tissue, 
        factor, 
        design, 
        include_factor_levels, 
        sample_metadata_full,
        SAMPLES_TO_EXCLUDE, 
        OUTPUT_DIR, 
        ALPHA, 
        LOG_CHANGE_CUTOFF,
        EXPRESSION_CUTOFF,
        EXPRESSION_CUTOFF_FRACTION
      )
      
    }
  }
}

# ======= Tests ========
# counts_files <- COUNTS_FILES
# tissue <- 'pit'
# factor <- 'month'
# design <- ~stage + month
# #include_factor_levels <- list('M', 'F', 'J')#list('F1', 'F2', 'F3', 'F4', 'F5')
# samples_to_exclude <- SAMPLES_TO_EXCLUDE
# output <- OUTPUT_DIR
# alpha <- ALPHA
# fold_change_cutoff <- LOG_CHANGE_CUTOFF
# expression_cutoff <- EXPRESSION_CUTOFF
# expression_cutoff_fraction <- EXPRESSION_CUTOFF_FRACTION
