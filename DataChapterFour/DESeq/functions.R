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


create_pca_plot <- function(dds, factor, tissue = '') {
  vsd <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup=c(factor), returnData=TRUE)
  colnames(pca_data)[colnames(pca_data) == factor] <- 'variable'
  pca_data$variable <- as.character(pca_data$variable)
  plt <- ggplot(pca_data, aes(PC1, PC2, color=variable, text=name)) + 
    geom_point() + scale_color_discrete(name = factor) + ggtitle(tissue)
  
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
    month = str_replace_all(metadata$Date, '/', '_')#str_replace_all(paste(metadata$Date, metadata$Sex), '/', '-')
  )
  
  return(df)
}
