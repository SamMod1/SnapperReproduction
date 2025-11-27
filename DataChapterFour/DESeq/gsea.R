setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
if (!require("pacman")) install.packages("pacman")
pacman::p_load(clusterProfiler, stringr, pathview, ggplot2)
if (!('enrichplot' %in% installed.packages())){BiocManager::install("enrichplot")}
library(enrichplot)
source("C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/GoEnrichment/go_enrichment.R") # Takes a little time to initialise but adds the function:


OUTPUT <- "DataChapterFour/DESeq/Results/SexGSEA"
DESEQ_RESULTS <- "DataChapterFour/DESeq/Results/"
FOLD_CHANGE_CUTOFF <- 1
ALPHA <- 0.05
SEED <- 1
MIN_NDEG <- 10
N_SHOWN <- 4


set.seed(SEED)

process_gene_list <- function(gene_list) {
  names(gene_list) <- str_replace_all(names(gene_list), regex('gene-'), '')
  
  multiple_homologues_bool <- grepl('_ChrAur', names(gene_list), fixed = TRUE)
  multiple_homologues <- names(gene_list)[multiple_homologues_bool]
  multiple_homologues <- str_replace_all(multiple_homologues, regex('_ChrAur[0-9]{0,9}'), '')
  multiple_homologues <- unique(multiple_homologues)
  
  new_gene_list <- c()
  
  for (homologue in multiple_homologues) {
    gene_mean <- c(mean(gene_list[grepl(homologue, names(gene_list), fixed = TRUE)]))
    names(gene_mean) <- c(homologue)
    new_gene_list <- append(new_gene_list, gene_mean)
  }
  
  new_gene_list <- new_gene_list[
    new_gene_list > FOLD_CHANGE_CUTOFF | new_gene_list < -FOLD_CHANGE_CUTOFF
  ]
  
  new_gene_list <- append(new_gene_list, gene_list[!multiple_homologues_bool])
  
  new_gene_list <- new_gene_list[!is.na(new_gene_list)]
  
  return(new_gene_list)
}


gene_list_from_file <- function(file, dir=DEG_DATA_DIR, sep=SEP, universe=FALSE) {
  path <- paste(dir, file, sep='/')
  df <- read.csv(path, sep=sep)
  if (universe) {return(df$gene)}
  df <- subset(
    df,
    df$padj <= ALPHA &
      (df$log2FoldChange > FOLD_CHANGE_CUTOFF | df$log2FoldChange < -FOLD_CHANGE_CUTOFF)
  )
  gene_list <- df$stat
  names(gene_list) <- df$gene
  
  return(gene_list)
}


modify_duplicates <- function(gene_list) {
  for (idx in which(duplicated(gene_list))) {
    gene_list[idx] <- gene_list[idx] + 0.000001
  }
  if (any(duplicated(gene_list))) {
    return(modify_duplicates(gene_list))
  } else(return(gene_list))
}


saved_wd <- getwd()
setwd(OUTPUT)

gene_mapping_df <- read.csv(GENE_MAPPING_FILE)[,c('X.tax_id', 'Symbol', 'description')]
gene_mapping_df <- gene_mapping_df[order(gene_mapping_df[,'Symbol']),]
gene_mapping_df <- gene_mapping_df[!duplicated(gene_mapping_df$Symbol),]
gene_mapping_df <- subset(gene_mapping_df, !substr(gene_mapping_df$description, 1, 15) == 'uncharacterized')
gene_mapping_df <- subset(gene_mapping_df, !gene_mapping_df$Symbol==gene_mapping_df$description)
gene_mapping <- gene_mapping_df$X.tax_id
names(gene_mapping) <- gene_mapping_df$Symbol

for (file in list.files(DESEQ_RESULTS)) {
  
  if (!grepl('.tsv', file, fixed = TRUE)){next}
  
  print(file)
  
  gene_list <- gene_list_from_file(file, DESEQ_RESULTS, sep='\t')
  go_list <- gene_list
  gene_list <- process_gene_list(gene_list)
  gene_list <- gene_list[names(gene_list)[names(gene_list) %in% names(gene_mapping)]]
  
  gene_ids <- unname(gene_mapping[names(gene_list)])
  names(gene_list) <- gene_ids
  
  if (any(duplicated(gene_list))) {gene_list <- modify_duplicates(gene_list)}
  if (any(duplicated(names(gene_list)))) {gene_list <- gene_list[!duplicated(names(gene_list))]}
  
  gene_list <- gene_list[order(gene_list, decreasing=TRUE)]
  
  if (length(gene_list) >= MIN_NDEG) {
    go_list <- go_list[order(go_list, decreasing=TRUE)]
    
    go <- gsea(
      go_list, 
      'GO', 
      title = file, 
      return_plot = FALSE, 
      alpha = 0.05,
      n_shown = 4
    )
    
    gse_kegg <- gseKEGG(
      gene_list,
      organism = 'dre',
      keyType = 'ncbi-geneid',
      pvalueCutoff = 0.05,
      verbose = FALSE,
      seed = SEED,
      #universe=gene_list_from_file(file, DEG_DATA_DIR, universe=TRUE)
    )
    
    if (!file.exists(substr(file, 1, nchar(file)-4))) {
      dir.create(substr(file, 1, nchar(file)-4))
    }
    
    setwd(substr(file, 1, nchar(file)-4))
    
    if (!is.null(go)) {
      try(dev.off(), silent=TRUE) 
      plt <- dotplot(go, showCategory=N_SHOWN, split='.sign', font.size = 10) + facet_grid(.~.sign)
      ggsave('GO.png', plot = plt)
      write.csv(go@result, 'GO.csv', row.names=FALSE)
    }
    
    gse_kegg@result <- subset(gse_kegg@result, gse_kegg@result$p.adjust < ALPHA)
    
    if (nrow(gse_kegg@result) > 0) {
      fig <- dotplot(gse_kegg, showCategory=N_SHOWN, split='.sign') + facet_grid(.~.sign)
      ggsave('KEGG_dotplot.png', plot=fig)
      write.csv(gse_kegg@result, 'KEGG.csv', row.names=FALSE)
    }
    
    if (
      nrow(gse_kegg@result) > 0
    ) {
      for (idx in 1:nrow(gse_kegg@result)) {
        row <- gse_kegg[idx,]
        
        print(row$ID)
        print(row$setSize)
        
        try(path <- pathview(
          gene.data = gene_list,
          pathway.id = print(row$ID),
          species = "dre",
          limit = list(gene=5, cpd=5),
          out.suffix = substr(file, 1, nchar(file) - 4)
        ), silent = TRUE)
      }
    }
    setwd('..')
  }
}

for (directory in list.files()) {
  if (length(list.files(directory)) == 0) {unlink(directory, force=TRUE, recursive = TRUE)}
}

setwd(saved_wd)
