setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # Only works when not running from source in RStudio
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggvenn, gplots, stringr, ggplot2, clusterProfiler, DBI, tools, UpSetR, pathview)
source("C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/GoEnrichment/go_enrichment.R")


ALPHA <- 0.05
EXPRESSION_CUTOFF <- 1
TISSUES <- c('gonad', 'pit', 'mid', 'brain')
DATA <- 'C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq2/Output - Sex V Month Wald/Data'
CONTRASTS <- list(c('F', 'J', 'M'), c('F1', 'F2', 'F3', 'F4', 'F5'), c('M2', 'M3', 'M4'))
#CONTRASTS <- list(c('17-11-2020 F', '09-06-2020 F', '06-10-2020 F'), c('17-11-2020 M', '09-06-2020 M', '06-10-2020 M'))
OUTPUT_VENN <- 'temp/UpSet'
OUTPUT_KEGG <- 'temp/OverRepresentation'
OUTPUT_DEGS <- "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/DESeq2/temp/DEGs"
SEP <- '\t'
GENE_MAPPING_FILE <- "C:/Users/cfnsjm/Local Doccuments/PhD/Programming/R/GoEnrichment/joined_gene_info.csv"
DR_DB <- 'C:/Users/cfnsjm/AppData/Local/Programs/R/R-4.4.2/library/org.Dr.eg.db/extdata/org.Dr.eg.sqlite'
PLOT_WIDTH <- 1080 
PLOT_HEIGHT <- 820


process_gene_list_kegg <- function(gene_list, gene_mapping) {
  genes <- sub('gene-', '', gene_list)
  genes <- sub(regex('_ChrAur[0-9]*'), '', genes)
  gene_ids <- unique(unname(gene_mapping[genes]))
  gene_ids <- gene_ids[!is.na(gene_ids)]
  
  return(gene_ids)
}


dr_db <- dbConnect(RSQLite::SQLite(), DR_DB)
kegg_ids <- dbGetQuery(dr_db, 'SELECT * FROM kegg')
ncbi_ids <- dbGetQuery(dr_db, 'SELECT * FROM genes')
dbDisconnect(dr_db)

gene_mapping_df <- read.csv(GENE_MAPPING_FILE)[,c('X.tax_id', 'Symbol', 'description')]
gene_mapping_df <- gene_mapping_df[order(gene_mapping_df[,'Symbol']),]
gene_mapping_df <- gene_mapping_df[!duplicated(gene_mapping_df$Symbol),]
gene_mapping_df <- subset(gene_mapping_df, !substr(gene_mapping_df$description, 1, 15) == 'uncharacterized')
gene_mapping_df <- subset(gene_mapping_df, !gene_mapping_df$Symbol==gene_mapping_df$description)
gene_mapping <- gene_mapping_df$X.tax_id
names(gene_mapping) <- gene_mapping_df$Symbol

males_all <- list()
females_all <- list()
all_ <- list()

if (!file.exists(OUTPUT_KEGG)) {dir.create(OUTPUT_KEGG)}
if (!file.exists(OUTPUT_VENN)) {dir.create(OUTPUT_VENN)}

files <- list.files(DATA)

origional_dir <- getwd()
tmp <- paste(origional_dir, 'TMP', sep='/')
dir.create(tmp)
setwd(OUTPUT_KEGG)
output_kegg <- getwd()
setwd(origional_dir)
setwd(OUTPUT_VENN)
output_venn <- getwd()
setwd(origional_dir)

for (tissue in TISSUES) {
  data <- list()
  tissue_files <- files[grepl(tissue, files)]
  gene_universe <- character()

  for (file in tissue_files) {
    if (!file_ext(file) == 'tsv' | grepl(',', file)) {next()}
    
    df <- read.csv(paste(DATA, file, sep='/'), sep=SEP)
    
    gene_universe <- append(gene_universe, df$gene)
    
    df <- subset(df, df$padj <= ALPHA)
    
    genes_up <- subset(df, df$log2FoldChange >= EXPRESSION_CUTOFF)$gene
    genes_down <- subset(df, df$log2FoldChange <= -EXPRESSION_CUTOFF)$gene
    
    if (grepl(regex('[A-Z][0-9]'), file)) {
      up <- substring(str_extract(file, regex('[A-Z][0-9] vs')), 1, 2)
      down <- substring(str_extract(file, regex('vs [A-Z][0-9]')), 4, 5)
    } else {
      up <- substring(str_extract(file, regex('[A-Z] vs')), 1, 1)
      down <- substring(str_extract(file, regex('vs [A-Z]')), 4, 4)
    }
    
    if (is.null(data[up][[1]])) {data[[up]][[1]] <- character()}
    if (is.null(data[down][[1]])) {data[down][[1]] <- character()}
    if (!is.null(genes_up)) {data[up][[1]] <- unique(append(data[up][[1]], genes_up))}
    if (!is.null(genes_down)) {data[down][[1]] <- unique(append(data[down][[1]], genes_down))}
  }
  
  gene_universe_kegg <- process_gene_list_kegg(gene_universe, gene_mapping)
  
  data <- data[!is.na(names(data))]
  
  for (condition in names(data)) {
    gene_ids_kegg <- process_gene_list_kegg(data[[condition]], gene_mapping)
    
    degs_filename <- paste(OUTPUT_DEGS, paste(condition, 'txt', sep='.'), sep='/')
    write(unlist(data[[condition]]), degs_filename, sep='\n')
    
    enriched_go <- enrichment(
      unlist(data[[condition]]),
      'GO',
      gene_universe = gene_universe,
      title = condition,
      return_plot = TRUE
    )
    
    enriched_kegg <- enrichKEGG(
      gene = as.character(gene_ids_kegg),
      organism = 'dre',
      keyType = 'kegg',
      pvalueCutoff = ALPHA,
      universe = as.character(gene_universe_kegg)
    )

    dir_name <- output_kegg %>% paste(
      tissue, sep='/') %>% paste(
        condition)
    
    if (!is.null(enriched_go) | !is.null(enriched_kegg)) {
      if (!dir.exists(dir_name)) {dir.create(dir_name)}
    }
    
    if (!is.null(enriched_go)) {try(dev.off(), silent=TRUE); ggsave(paste(dir_name, 'GO.png', sep='/'), enriched_go)}

    if (FALSE) {#(!is.null(enriched_kegg)) {
      if (nrow(enriched_kegg) > 0) {
        plt <- dotplot(enriched_kegg, showCategory=10, title=paste(tissue, condition))
        ggsave(paste(dir_name, 'barplot.png', sep='/'), plt)

        setwd(tmp)

        gene_list <- rep(5, length(gene_ids_kegg))
        names(gene_list) <- gene_ids_kegg

        results <- subset(enriched_kegg@result, enriched_kegg@result$p.adjust <= ALPHA)

        for (idx in 1:nrow(results)) {
          row <- results[idx,]

          try(path <- pathview(
            gene.data = gene_list,
            pathway.id = row$ID,
            species = "dre",
            limit = list(gene=5, cpd=5),
            out.suffix = paste(tissue, condition)
          ), silent = TRUE)
          out_file <- paste(row$ID, paste(tissue, condition), sep='.') %>%
            paste('png', sep='.')
          file.rename(paste(tmp, out_file, sep='/'), paste(dir_name, out_file, sep='/'))
        }
        while (!is.null(dev.list()))  dev.off()
      }
    }
  }
  
  for (contrast in CONTRASTS) {
    venn_data <- data[names(data) %in% contrast]
    # plt <- ggvenn(venn_data)
    plt_name <- output_venn %>% paste(
      tissue, sep='/'
      ) %>% paste(
      toString(contrast), sep=' '
      ) %>% paste(
      'png', sep='.'
      )
    #png(file = plt_name, width = PLOT_WIDTH, height = PLOT_HEIGHT)
    #print(upset(
    #  fromList(data[names(data) %in% contrast]), 
    #  order.by = 'freq', 
    #  show.numbers = FALSE, 
    #  text.scale = 2, 
    #))
    #dev.off()
    #ggsave(plt_name)
}
  
  tissue_all <- character()
  tissue_males <- character()
  tissue_females <- character()
  
  for (condition in names(data)) {
    if (substring(condition, 1, 1) == 'M') {
      tissue_males <- append(tissue_males, data[[condition]])
    } else if (substring(condition, 1, 1) == 'F') {
      print(condition)
      tissue_females <- append(tissue_females, data[[condition]])
      print(length(tissue_females))
    }
    tissue_all <- append(tissue_all, data[[condition]])
  }
  all_[[tissue]] <- unique(tissue_all)
  males_all[[tissue]] <- unique(tissue_males)
  females_all[[tissue]] <- unique(tissue_females)
}

# plt_name <- paste(output_venn, 'All tissues.png', sep='/')
# png(filename = plt_name, width = PLOT_WIDTH, height = PLOT_HEIGHT)
# print(upset(fromList(all_), order.by = 'freq'))
# dev.off()#ggsave(plt_name, plt)
# 
plt_name <- paste(output_venn, 'All tissues Female.png', sep='/')
png(filename = plt_name, width = PLOT_WIDTH, height = PLOT_HEIGHT)
#upset(fromList(females_all), order.by = 'freq')
ggvenn(females_all, text_size = 14, show_percentage = FALSE, stroke_size = 2)
dev.off()#ggsave(plt_name, plt)
# 
plt_name <- paste(output_venn, 'All tissues Male.png', sep='/')
png(filename = plt_name, width = PLOT_WIDTH, height = PLOT_HEIGHT)
#upset(fromList(males_all), order.by = 'freq')
ggvenn(males_all, text_size = 14, show_percentage = FALSE, stroke_size = 2)
dev.off()#ggsave(plt_name, plt)

setwd(origional_dir)
unlink(tmp, recursive=TRUE)
