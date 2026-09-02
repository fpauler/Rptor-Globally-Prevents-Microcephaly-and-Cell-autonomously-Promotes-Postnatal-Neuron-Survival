# =============================================
# sixth step of analysis
# script to perform a STRING analysis of DEPs
# =============================================

library (openxlsx)
library (STRINGdb)
library (igraph)
library (RCy3)
library (pheatmap)
library (tidyverse)
library (ggupset)

# =======================
# set up I/O folders
# =======================
base <- "./bulk_RNASeq/biorxiv/"

# ==========================================================
# read the data
# Note: I focus this analysis on DEPs with GO associations
# ==========================================================
all_dep_df <- data.frame()
eGO <- data.frame()

for (age in c("E13", "E16", "P0")) {
  # collect all genes via GO terms
  # for all developmental time points (ages) and only for Rptor-MADM (het)
  tmp_eGO <- read.xlsx( paste(base, "Supplement/DEP_GO_enrichment.xlsx", sep=""), sheet = paste(age, "het", sep="_") )
  
  tmp_eGO$age <- age
  eGO <- rbind( eGO, tmp_eGO)
  
  all_dep_df <- rbind( all_dep_df, data.frame( age = age, symbol = unique( unlist( sapply( tmp_eGO$geneID, function (x) strsplit(x, "/")[[1]] ) ))))
}

table(table(all_dep_df$symbol)) 
#   1   2   3 
# 314  52   1 

# I only use dev. time specific genes!
# in other words: 53 genes are excluded in this analysis
# this is the gene group I will use for constructing the network
all_dep_df <- all_dep_df[ which( !(duplicated(all_dep_df$symbol) | duplicated(all_dep_df$symbol, fromLast = TRUE))), ]

# write out DEP list that goes into STRING network analysis
wb <- createWorkbook()
tmp_df <- all_dep_df
sheetName <- "DEP for STRING"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, tmp_df )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( tmp_df ) )
setColWidths( wb, sheetName, cols = 1:ncol( tmp_df ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/DEP_STRING_raw.xlsx", sep=""), overwrite = T) 


# =======================================================
# create a simple GO term grouping for simpler plotting
# this part creates a map of genes to GO groups
# =======================================================

eGO$group <- ""
for (group in c("axon", "cell_cycle", "mitoch", "TOR", "apop", "reactive_oxygen", "metabol", "cytokine")) {
  eGO$group[ which(grepl(pattern = group, x = eGO$Description)) ]  <- group
}

# two terms are not used in this analysis
eGO$Description[ which(eGO$group == "") ]
# [1] "regulation_of_necroptotic_process" "necroptotic_process"

eGO <-  eGO[ which(!eGO$group == ""), ]
gene2group <- sapply( unique(eGO$group), function (group){
  tmp_df <- eGO[ which(eGO$group == group),]
  return ( unique( unlist( sapply( tmp_df$geneID, function (x) strsplit(x, "/")[[1]] ) )) )
})

# extract unique genes
ambiguous_genes <- names( which( table(unlist(gene2group)) > 1 ) )

# lookup vector
gene_vec <- c()
name_vec <- c()

for (name in names(gene2group)) {
  tmp_gene_list <- setdiff(gene2group[[name]], ambiguous_genes)
  gene_vec <- c(gene_vec, tmp_gene_list )
  name_vec <- c(name_vec, rep(name, length(tmp_gene_list)) )
}

gene2group_lookup <- name_vec
names( gene2group_lookup ) <- gene_vec

# add multiple matchers
gene2group_lookup <- c( name_vec, rep( "multiple", length(ambiguous_genes)) )
names( gene2group_lookup ) <- c( gene_vec, ambiguous_genes )

# ==========================================================================
# STRING analysis
# Note: Cytoscape network used for final figure presentation available as:
#       DEG_DEP_STRING.cys on GitHub
# ==========================================================================

# STRINGdb analysis  -  downloaded 2.5.2025
string_db <- STRINGdb$new( version="12", species=10090, score_threshold=400, input_directory=paste(base, "other_data", sep=""))

# map symbols to STRING and get interactions
genes_4_nw <- string_db$map( all_dep_df, removeUnmappedRows = TRUE, my_data_frame_id_col_names = "symbol")
complete_network <- string_db$get_interactions( genes_4_nw$STRING_id )

pdf(file = paste(base, "/plots/QC/STRING_DEP.pdf", sep=""))
string_db$plot_network( genes_4_nw$STRING_id )
dev.off()

g <- graph_from_data_frame(complete_network, directed=F)

#remove multiple edges
g_simple <- delete_edges(g, E(g)[which_multiple(g)])

g_simple_clusters <- components(g_simple)

#get the largest connected cluster in the network
idx <- which(g_simple_clusters$csize == max(g_simple_clusters$csize))
connected_genes <- names(which(g_simple_clusters$membership == idx))

length ( connected_genes )
# [1] 285

largest_connected_nw <- string_db$get_interactions( connected_genes )
g <- graph_from_data_frame(largest_connected_nw, directed=F)

#remove multiple edges
g_simple <- delete_edges(g, E(g)[which_multiple(g)])

V(g_simple)$name <- sapply( V(g_simple)$name, function (x) if ( x == '') {return('')} else {genes_4_nw$symbol[which(genes_4_nw$STRING_id == x)]})
V(g_simple)$group <- sapply( V(g_simple)$name, function (x) as.character( gene2group_lookup[x] ) )
V(g_simple)$group <- sapply( V(g_simple)$group, function (x) if (any(is.na(x))) { return("other")} else {return (x)} )
V(g_simple)$age <- sapply( V(g_simple)$name, function (x) { all_dep_df$age[which(all_dep_df$symbol == x)]})

c("axon", "cell_cycle", "mitoch", "TOR", "apop", "reactive_oxygen", "metabol", "cytokine")

group2color = c(  "axon" = "#0000FF", 
                  "mitoch" = "#FF8C00", 
                  "cytokine" = "#00FF00", 
                  "apop" = "#C1CDC1", 
                  "reactive_oxygen" = "#9400D3", 
                  "metabol" = "#00CED1",
                  "multiple" = "#CDCD00")

V(g_simple)$color <- group2color[V(g_simple)$group]
V(g_simple)$color[ which(is.na(V(g_simple)$color))] <- "#FFFAF0" 
V(g_simple)$font_size <- 0

V(g_simple)$font_size [ which( V(g_simple)$name %in% c(names(gene2group_lookup), ambiguous_genes )) ] <- 12

createNetworkFromIgraph( g_simple ,"DEP network 1")

# upset plot to make a point about functional connections

group_mat <- sapply( V(g_simple)$name, function (gene){
  sapply( c("TOR", "cell_cycle", "apop", "axon", "reactive_oxygen", "metabol", "mitoch", "cytokine") , function (group){
    ifelse ( gene %in% gene2group[[group]], 1, 0 )
  })
})

# limit complexity of the plot - focus on 4 GO term groups
connect_gene_df <- reshape2::melt( gene2group[ c("mitoch", "apop", "axon", "reactive_oxygen") ] )
connect_gene_df %>% group_by(value) %>% summarize(pathway = list(L1), n=nchar(pathway) ) -> connect_gene_tibble
connect_gene_tibble$comb <- sapply(connect_gene_tibble$pathway, function (x) paste(sort(x), collapse=","))

connect_gene_tibble %>%
  ggplot(aes(x = pathway)) + geom_bar() + scale_x_upset() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  theme_classic() + ylim(0,130) -> upset_plot

ggsave( filename = paste(base, "plots/Main/Fig3I_DEP_STRING_upset_plot.pdf", sep=""), plot = upset_plot )

# write out GO classification - complete list is given
wb <- createWorkbook()
tmp_df <- as.data.frame( connect_gene_tibble[,c("value", "comb")] )
sheetName <- "GO group association"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, tmp_df )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( tmp_df ) )
setColWidths( wb, sheetName, cols = 1:ncol( tmp_df ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/DEP_STRING_GO_classification.xlsx", sep=""), overwrite = T) 

sessionInfo()

# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.3.2/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8    
# [5] LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggupset_0.4.1    lubridate_1.9.3  forcats_1.0.0    stringr_1.5.1    dplyr_1.1.4      purrr_1.0.2      readr_2.1.5     
# [8] tidyr_1.3.1      tibble_3.2.1     ggplot2_3.5.0    tidyverse_2.0.0  pheatmap_1.0.12  RCy3_2.22.1      igraph_2.1.4    
# [15] STRINGdb_2.14.3  openxlsx_4.2.5.2
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1    IRdisplay_1.1       farver_2.1.1        blob_1.2.4          bitops_1.0-7        fastmap_1.1.1      
# [7] RCurl_1.98-1.14     XML_3.99-0.16.1     digest_0.6.34       base64url_1.4       timechange_0.3.0    lifecycle_1.0.4    
# [13] RSQLite_2.3.6       magrittr_2.0.3      compiler_4.3.2      rlang_1.1.3         tools_4.3.2         plotrix_3.8-4      
# [19] utf8_1.2.4          labeling_0.4.3      bit_4.0.5           curl_5.2.0          plyr_1.8.9          repr_1.1.7         
# [25] RColorBrewer_1.1-3  KernSmooth_2.23-22  pbdZMQ_0.3-14       withr_3.0.0         BiocGenerics_0.48.1 grid_4.3.2         
# [31] hash_2.2.6.3        stats4_4.3.2        fansi_1.0.6         caTools_1.18.2      colorspace_2.1-0    scales_1.3.0       
# [37] gtools_3.9.5        cli_3.6.2           crayon_1.5.2        chron_2.3-62        generics_0.1.3      rstudioapi_0.16.0  
# [43] reshape2_1.4.4      httr_1.4.7          tzdb_0.4.0          DBI_1.2.2           cachem_1.0.8        base64enc_0.1-3    
# [49] vctrs_0.6.5         jsonlite_1.8.8      hms_1.1.3           bit64_4.0.5         packrat_0.9.2       proto_1.0.0        
# [55] glue_1.7.0          RJSONIO_2.0.0       stringi_1.8.3       gtable_0.3.4        munsell_0.5.0       pillar_1.9.0       
# [61] htmltools_0.5.7     gplots_3.1.3.1      IRkernel_1.3.2      graph_1.80.0        R6_2.5.1            evaluate_0.23      
# [67] gsubfn_0.7          png_0.1-8           backports_1.4.1     memoise_2.0.1       Rcpp_1.0.12         zip_2.3.1          
# [73] uuid_1.2-0          sqldf_0.4-11        fs_1.6.3            pkgconfig_2.0.3
