setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggvenn, gplots, stringr, ggplot2, clusterProfiler, DBI, tools, UpSetR, pathview)
source("C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/GoEnrichment/go_enrichment.R")


FOLD_CHANGE_CUTOFF <- 1
ALPHA <- 0.05
OUTPUT <- "DataChapterTwo/DESeq/Results/"
DESEQ_RESULTS <- "DataChapterTwo/DESeq/Results/F" #Females
# DESEQ_RESULTS <- "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq_final/Results/M" #Males


genes_from_file <- function(file, dir=DESEQ_RESULTS, sep='\t', universe=FALSE) {
  path <- paste(dir, file, sep='/')
  df <- read.csv(path, sep=sep)
  if (universe) {return(df$gene)}
  df <- subset(
    df,
    df$padj <= ALPHA &
      (df$log2FoldChange > FOLD_CHANGE_CUTOFF | df$log2FoldChange < -FOLD_CHANGE_CUTOFF)
  )
  gene_list_up <- df$gene[df$log2FoldChange > 0]
  gene_list_down <- df$gene[df$log2FoldChange < 0]
  
  return(list(gene_list_up, gene_list_down))
}


process_gene_list <- function(gene_list) {
  gene_list <- str_replace_all(gene_list, 'gene-', '')
  #gene_list <- str_replace_all(gene_list, regex('_ChrAur[0-9]{0,9}'), '')
  gene_list <- unique(gene_list)
  gene_list <- gene_list[!is.na(gene_list)]
  return(gene_list)
}


stage_list <- list()

for (file in list.files(DESEQ_RESULTS)) {
  
  if (!grepl('.tsv', file, fixed = TRUE)){next}
  
  contrast_1 <- str_split(file, ' vs ')[[1]][1]
  contrast_2 <- substr(str_split(file, ' vs ')[[1]][2], 1, 2)
  
  genes <- genes_from_file(file)
  genes_1 <- process_gene_list(genes[[1]])
  genes_2 <- process_gene_list(genes[[2]])
  
  stage_list[[contrast_1]] <- append(stage_list[[contrast_1]], genes_1)
  stage_list[[contrast_2]] <- append(stage_list[[contrast_2]], genes_2)
}

for (stage in names(stage_list)) {stage_list[[stage]] <- unique(stage_list[[stage]])}
for (stage in names(stage_list)) {print(stage);print(length(stage_list[[stage]]))}
#plt_name <- paste(output_venn, 'All tissues Female.png', sep='/')
#png(filename = plt_name, width = PLOT_WIDTH, height = PLOT_HEIGHT)
upset(fromList(stage_list), order.by='freq', text.scale=1.5)#, empty.intersections=TRUE)

#dev.off()#ggsave(plt_name, plt)
