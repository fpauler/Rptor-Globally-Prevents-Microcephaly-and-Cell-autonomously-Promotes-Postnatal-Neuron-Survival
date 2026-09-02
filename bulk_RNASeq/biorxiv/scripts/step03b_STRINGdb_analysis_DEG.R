# ============================
# fourth step of analysis
# this script does:
# plot basic DEG stats
# GO term enrichment
# STRING analysis for DEGs
# ==============================

library (openxlsx)
library (STRINGdb)
library (igraph)
library (RCy3)
library (clusterProfiler)
library (ggplot2)
library (ggVennDiagram)
library (cowplot)
library (tidyverse)
library (ggupset)

set.seed(2401)

# =======================
# set up I/O folders
# =======================
base <- "./bulk_RNASeq/biorxiv/"

# load filtered GO term df
GO_GENE_DF <- readRDS( file = paste(base, "/other_data/GO_Term_df.RDS", sep="") )
GO_GENE_DF <- GO_GENE_DF[,c(2,1)]

# basic DEG stats
eGO <- list()
deg_list <- list()
deg_count <- list()

# ======================================================================
# read the DEG analysis, filter DEGs and do GO term enrichment analysis
# ======================================================================

for (age in c("E13", "E16", "P0")) {
  for ( group in c("cKO", "HET") ) {
    deg_df <- read.xlsx( paste(base, "Supplement/DEG_raw_DEG_stats.xlsx", sep=""), sheet = paste( age, group, "WT", sep="_" ))
    
    deg_list[[ paste( age, group, "all", sep="_") ]] <- deg_df$symbol[ which(deg_df$padj < 0.05) ]
    deg_list[[ paste( age, group, "up", sep="_") ]] <- deg_df$symbol[ which(deg_df$padj < 0.05 & deg_df$log2FoldChange > 0) ]
    deg_list[[ paste( age, group, "down", sep="_") ]] <- deg_df$symbol[ which(deg_df$padj < 0.05 & deg_df$log2FoldChange < 0) ]
  
    eGO[[ paste( age, group, sep="_") ]] <- NA
    eGO[[ paste( age, group, sep="_") ]] <-  enricher ( gene = deg_list[[ paste( age, group, "all", sep="_") ]], 
                                                           TERM2GENE = GO_GENE_DF, pvalueCutoff = 0.9, qvalueCutoff = 0.9,
                                                           minGSSize = 10, 
                                                           maxGSSize = 2000
    )
  }
}


# write out enrichment analysis
# this file contains all GO enrichment output
wb <- createWorkbook()
for (name in names(eGO)) {
  
  sheetName <- name
  addWorksheet(wb, sheetName)
  tmp_df <- eGO[[ sheetName ]]
  
  writeData( wb, sheetName, tmp_df )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol( tmp_df ) )
  setColWidths( wb, sheetName, cols = 1:ncol( tmp_df ), widths="auto" )
}
saveWorkbook(wb, file = paste(base, "Supplement/DEG_GO_enrichment.xlsx", sep=""), overwrite = T) 

# DEG count figure

deg2plot <- reshape2::melt( sapply(deg_list, length ) )
deg2plot$age <- sapply( rownames(deg2plot), function (x) strsplit(x,"_")[[1]][1])
deg2plot$genotype <- sapply( rownames(deg2plot), function (x) strsplit(x,"_")[[1]][2])
deg2plot$dir <- sapply( rownames(deg2plot), function (x) strsplit(x,"_")[[1]][3])
deg2plot$genotype <- factor( deg2plot$genotype, levels = c("HET",  "cKO"))

ggplot(deg2plot[which(!deg2plot$dir == "all"),], aes(x=genotype, y=value, fill = dir)) + geom_bar(stat="identity") + 
  facet_wrap(~age) + theme_classic() + scale_fill_manual( values = c("up" = "goldenrod", "down" = "cornflowerblue")) + 
  ggtitle("#DEGs FDR 5%, relative to Control-MADM") + xlab("") + ylab("")

ggsave  ( filename = paste(base, "/plots/Main/Fig3B_DEG_number_bar_plot.pdf", sep=""), width=4, height=4)

# plot overlap cKO / HET for each age

plot_list <- list()
for (age in c("E13", "E16", "P0")) {
  
  tmp_deg_list <- list(c(deg_list[[ paste(age, "cKO_all", sep="_") ]]), c(deg_list[[ paste(age, "HET_all", sep="_") ]]) )
  names(tmp_deg_list) <- c(paste(age, "cKO", sep="_"), paste(age, "het", sep="_") )
  plot_list[[ age ]] <- ggVennDiagram( x = tmp_deg_list, label_geom = "label", label = "count", label_alpha = 0 ) + 
    scale_fill_gradient(low="white",high = "white") +
    theme(legend.position = "NONE") + coord_cartesian()
  
}

venn_plot <- plot_grid( plotlist = plot_list, ncol = 1, greedy = F )
ggsave  ( filename = paste(base, "/plots/Supplement/FigS3C_DEG_venn_plot_comp.pdf", sep=""), plot = venn_plot, width=6, height=20)

# plot overlap for ages separately for cKO / HET

plot_list <- list()
for (comp in c("cKO", "HET")) {
  
  tmp_deg_list <- list(c(deg_list[[ paste("E13", comp, "all", sep="_") ]]), 
                       c(deg_list[[ paste("E16", comp, "all", sep="_") ]]), 
                       c(deg_list[[ paste("P0", comp, "all", sep="_") ]]) )
  
  names(tmp_deg_list) <- c("E13", "E16", "P0")
  
  plot_list[[ comp ]] <- ggVennDiagram( x = tmp_deg_list, label = "count", label_geom = "label" ) + 
    scale_fill_gradient(low="white",high = "white") + theme(legend.position = "NONE") + coord_cartesian() +
    ggtitle( comp )
  
}

venn_plot <- plot_grid( plotlist = plot_list, ncol = 1, greedy = F )
ggsave  ( filename = paste(base, "/plots/Main/Fig3D_DEG_venn_plot_ages.pdf", sep=""), plot = venn_plot, width=3, height=5)

# write out number of DEGs - raw data for Fig3B
wb <- createWorkbook()

sheetName <- "DEG number"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, deg2plot )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( deg2plot ) )
setColWidths( wb, sheetName, cols = 1:ncol( deg2plot ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/Fig3B_DEG_numbers.xlsx", sep=""), overwrite = T) 

# ==================
# GO term analysis
# ==================

GO_vec <- c()
for (name in names(eGO)) {
  # extract significantly enriched terms and combine data in a way I can plot a heatmap later
  # NOTE: use pvalue here to get inclusive picture
  GO_vec <- c( GO_vec, eGO[[name]]@result$ID[which(eGO[[name]]@result$pvalue < 0.05)] )
}

# group terms based on Description
df <- data.frame( Description = unique(GO_vec), group ="")

for (group in c("axon", "cell_cycle", "mitoch", "TOR", "apop", "reactive_oxygen", "metabol", "cytokine")) {
  df$group[ which(grepl(pattern = group, x = df$Description)) ]  <- group
}
rownames(df) <- df$Description

# write out GO classification determined above
wb <- createWorkbook()

sheetName <- "GO groups"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, df )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( df ) )
setColWidths( wb, sheetName, cols = 1:ncol( df ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/DEG_GO_classification.xlsx", sep=""), overwrite = T) 

# group GO terms and plot max score per group
GOgroup_df2plot <- data.frame()
for (i in 1:length(eGO)) {
  
  id_parts <- strsplit( names(eGO)[i], "_")[[1]]
  
  tmp <- eGO[[i]]@result[which( eGO[[i]]@result$Description %in% df$Description ),]
  
  tmp$pvalue[ which( is.na(tmp$pvalue)) ] <- 1
  tmp$score <- log10(tmp$pvalue) * -1

  for (group in unique(df$group) ) {
    tmp_idx <- which( tmp$ID %in% df[which(df$group == group), "Description"] )
    if (length(tmp_idx) > 0) {
      max_score <- max( tmp$score[ tmp_idx ] )
    } else {
      max_score <- 0
    }
    
    GOgroup_df2plot <- rbind( GOgroup_df2plot, data.frame( GO_group = group, value = max_score, age = id_parts[1], genotype = id_parts[2]) )
  }

}

GOgroup_df2plot$value[which(GOgroup_df2plot$value > 2.5)] <- 2.5 
GOgroup_df2plot$GO_group <- factor( GOgroup_df2plot$GO_group, levels = rev( c("cell_cycle", "cytokine", "TOR", "reactive_oxygen", "metabol", "mitoch", "axon", "apop")))
GOgroup_df2plot$genotype <- factor( GOgroup_df2plot$genotype, levels = c("HET", "cKO"))

ggplot(GOgroup_df2plot, aes(x=genotype, y=GO_group, fill=value)) + geom_tile() + 
  # scale_fill_gradient2(low="blue", mid="white", high = "gold", midpoint = 1, limits = c(0,2.5)) + 
  scale_fill_gradient(low="white", high = "black", limits = c(0,2.5)) + 
  facet_wrap(~age) + theme_classic() + ggtitle("transcriptomic GO term analysis")
ggsave  ( filename = paste(base, "/plots/Supplement/FigS3D_DEG_GO_score_plot.pdf", sep=""), width=4, height=4)

# ==========================================================================
# STRING analysis
# Note: Cytoscape network used for final figure presentation available as:
#       DEG_DEP_STRING.cys on GitHub
# ==========================================================================

# focus on HET analysis - collect again all DEGs to test for connection throughout development
all_deg_df <- data.frame()
for (age in c("E13", "E16", "P0")) {
  deg_df <- read.xlsx( paste(base, "Supplement/DEG_raw_DEG_stats.xlsx", sep=""),sheet = paste( age, "HET_WT", sep="_" ))
  tmp_df <- deg_df[ which(deg_df$padj < 0.05), ]
  tmp_df$age <- age
  all_deg_df <- rbind( all_deg_df, tmp_df )
}

# confirm overlap analysis
table(table(all_deg_df$symbol))

#   1   2 
# 212   3 

# associate one/more functions with genes using GO  
all_deg_df <- all_deg_df[ which( !(duplicated(all_deg_df$symbol) | duplicated(all_deg_df$symbol, fromLast = TRUE))), ]

# staring number of DEGs
nrow(all_deg_df)
# [1] 212

# more classification than enrichment
eGO <- enricher ( gene = all_deg_df$symbol, TERM2GENE = GO_GENE_DF, pvalueCutoff = 0.9, qvalueCutoff = 0.9,
                  minGSSize = 10, 
                  maxGSSize = 2000
)

eGO <- eGO@result
eGO$group <- ""
for (group in c("axon", "cell_cycle", "mitoch", "TOR", "apop", "reactive_oxygen", "metabol", "cytokine")) {
  eGO$group[ which(grepl(pattern = group, x = eGO$Description)) ]  <- group
}

gene2group <- sapply( unique(eGO$group), function (group){
  tmp_df <- eGO[ which(eGO$group == group),]
  return ( unique( unlist( sapply( tmp_df$geneID, function (x) strsplit(x, "/")[[1]] ) )) )
})

# How many genes can be associated with a GO group?
length( unique( unlist(gene2group) ) )
# [1] 89

# extract unique genes
ambiguous_genes <- names( which( table(unlist(gene2group)) > 1 ) )
length(ambiguous_genes)
# [1] 34

# lookup vector
gene_vec <- c()
name_vec <- c()

for (name in names(gene2group)) {
  tmp_gene_list <- setdiff(gene2group[[name]], ambiguous_genes)
  gene_vec <- c(gene_vec, tmp_gene_list )
  name_vec <- c(name_vec, rep(name, length(tmp_gene_list)) )
}

# add multiple matchers
gene2group_lookup <- c( name_vec, rep( "multiple", length(ambiguous_genes)) )
names( gene2group_lookup ) <- c( gene_vec, ambiguous_genes )


# STRINGdb analysis  -  downloaded 2.5.2025
string_db <- STRINGdb$new( version="12", species=10090, score_threshold=400, input_directory=paste(base, "other_data", sep=""))

# map symbols to STRING and get interactions
genes_4_nw <- string_db$map( all_deg_df, removeUnmappedRows = TRUE, my_data_frame_id_col_names = "symbol")
# Warning:  we couldn't map to STRING 5% of your identifiers> 

complete_network <- string_db$get_interactions( genes_4_nw$STRING_id )

pdf(file = paste(base, "/plots/QC/STRING_DEG.pdf", sep=""))
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
# [1] 79

largest_connected_nw <- string_db$get_interactions( connected_genes )
g <- graph_from_data_frame(largest_connected_nw, directed=F)

#remove multiple edges
g_simple <- delete_edges(g, E(g)[which_multiple(g)])

V(g_simple)$name <- sapply( V(g_simple)$name, function (x) if ( x == '') {return('')} else {genes_4_nw$symbol[which(genes_4_nw$STRING_id == x)]})
V(g_simple)$group <- sapply( V(g_simple)$name, function (x) as.character( gene2group_lookup[x] ) )
V(g_simple)$group <- sapply( V(g_simple)$group, function (x) if (any(is.na(x))) { return("other")} else {return (x)} )
V(g_simple)$age <- sapply( V(g_simple)$name, function (x) { all_deg_df$age[which(all_deg_df$symbol == x)]})

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

# How many genes in the network:
length( V(g_simple)$group )
# [1] 79

# send to Cytoscape for making figure
createNetworkFromIgraph( g_simple ,"DEG network")

# prepare an overview plot showing GO connection between genes
connect_gene_df <- reshape2::melt( gene2group[ c("mitoch", "apop", "axon", "reactive_oxygen") ] )
connect_gene_df %>% group_by(value) %>% summarize(pathway = list(L1), n=nchar(pathway) ) -> connect_gene_tibble
connect_gene_tibble$comb <- sapply(connect_gene_tibble$pathway, function (x) paste(sort(x), collapse=","))

connect_gene_tibble %>%
  ggplot(aes(x = pathway)) + geom_bar() + scale_x_upset() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1) +
  theme_classic() + ylim(0,20) -> upset_plot
ggsave( filename = paste(base, "plots/Main/Fig3H_DEG_STRING_upset_plot.pdf", sep=""), plot = upset_plot )

# write out GO connection of genes
wb <- createWorkbook()
tmp_df <- as.data.frame( connect_gene_tibble[,c("value", "comb")] )
sheetName <- "GO group association"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, tmp_df )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( tmp_df ) )
setColWidths( wb, sheetName, cols = 1:ncol( tmp_df ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/DEG_STRING_GO_classification.xlsx", sep=""), overwrite = T) 

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
#   [1] ggupset_0.4.1          lubridate_1.9.3        forcats_1.0.0          stringr_1.5.1          dplyr_1.1.4           
# [6] purrr_1.0.2            readr_2.1.5            tidyr_1.3.1            tibble_3.2.1           tidyverse_2.0.0       
# [11] pheatmap_1.0.12        cowplot_1.1.3          ggVennDiagram_1.5.2    ggplot2_3.5.0          clusterProfiler_4.10.1
# [16] RCy3_2.22.1            igraph_2.1.4           STRINGdb_2.14.3        openxlsx_4.2.5.2      
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3         
# [5] farver_2.1.1            ragg_1.3.0              fs_1.6.3                zlibbioc_1.48.2        
# [9] vctrs_0.6.5             memoise_2.0.1           RCurl_1.98-1.14         ggtree_3.10.1          
# [13] base64enc_0.1-3         htmltools_0.5.7         curl_5.2.0              plotrix_3.8-4          
# [17] gridGraphics_0.5-1      KernSmooth_2.23-22      gsubfn_0.7              plyr_1.8.9             
# [21] cachem_1.0.8            uuid_1.2-0              lifecycle_1.0.4         pkgconfig_2.0.3        
# [25] gson_0.1.0              Matrix_1.6-5            R6_2.5.1                fastmap_1.1.1          
# [29] GenomeInfoDbData_1.2.11 digest_0.6.34           aplot_0.2.2             enrichplot_1.22.0      
# [33] colorspace_2.1-0        patchwork_1.2.0         AnnotationDbi_1.64.1    S4Vectors_0.40.2       
# [37] chron_2.3-62            textshaping_0.3.7       RSQLite_2.3.6           base64url_1.4          
# [41] labeling_0.4.3          timechange_0.3.0        fansi_1.0.6             RJSONIO_2.0.0          
# [45] httr_1.4.7              polyclip_1.10-6         compiler_4.3.2          bit64_4.0.5            
# [49] withr_3.0.0             backports_1.4.1         BiocParallel_1.36.0     viridis_0.6.5          
# [53] DBI_1.2.2               ggforce_0.4.2           gplots_3.1.3.1          MASS_7.3-60            
# [57] HDO.db_0.99.1           gtools_3.9.5            caTools_1.18.2          tools_4.3.2            
# [61] ape_5.8                 scatterpie_0.2.2        zip_2.3.1               glue_1.7.0             
# [65] nlme_3.1-163            GOSemSim_2.28.1         grid_4.3.2              shadowtext_0.1.3       
# [69] pbdZMQ_0.3-14           reshape2_1.4.4          fgsea_1.28.0            generics_0.1.3         
# [73] gtable_0.3.4            tzdb_0.4.0              hms_1.1.3               data.table_1.15.0      
# [77] tidygraph_1.3.1         utf8_1.2.4              XVector_0.42.0          BiocGenerics_0.48.1    
# [81] ggrepel_0.9.5           pillar_1.9.0            yulab.utils_0.1.4       IRdisplay_1.1          
# [85] splines_4.3.2           tweenr_2.0.3            treeio_1.26.0           lattice_0.21-9         
# [89] bit_4.0.5               tidyselect_1.2.1        GO.db_3.18.0            Biostrings_2.70.3      
# [93] gridExtra_2.3           IRanges_2.36.0          sqldf_0.4-11            stats4_4.3.2           
# [97] graphlayouts_1.1.1      Biobase_2.62.0          proto_1.0.0             stringi_1.8.3          
# [101] lazyeval_0.2.2          ggfun_0.1.4             evaluate_0.23           codetools_0.2-19       
# [105] ggraph_2.2.1            qvalue_2.34.0           hash_2.2.6.3            graph_1.80.0           
# [109] ggplotify_0.1.2         cli_3.6.2               IRkernel_1.3.2          systemfonts_1.0.6      
# [113] repr_1.1.7              munsell_0.5.0           Rcpp_1.0.12             GenomeInfoDb_1.38.8    
# [117] png_0.1-8               XML_3.99-0.16.1         parallel_4.3.2          blob_1.2.4             
# [121] DOSE_3.28.2             bitops_1.0-7            viridisLite_0.4.2       tidytree_0.4.6         
# [125] scales_1.3.0            packrat_0.9.2           crayon_1.5.2            rlang_1.1.3            
# [129] fastmatch_1.1-4         KEGGREST_1.42.0  
