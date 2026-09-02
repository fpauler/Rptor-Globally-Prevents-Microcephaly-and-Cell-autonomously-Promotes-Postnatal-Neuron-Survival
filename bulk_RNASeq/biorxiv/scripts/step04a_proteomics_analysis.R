# ====================================================================
# fifth step of analysis
# analysis of DEP
# differential abundance analysis from ISTA proteomics core facility
# ====================================================================

library(openxlsx)
library(ggplot2)
library(clusterProfiler)
library (cowplot)
library (ggVennDiagram)

set.seed(2401)

# =======================
# set up I/O folders
# =======================
base <- "./bulk_RNASeq/biorxiv/"

# =================================================================================
# define lookup list to where to find relevant datasets for each age and genotype
# =================================================================================
lookup_list <- list()
lookup_list[["E13"]] <- list()
lookup_list[["E13"]][["cKO"]] <- "T-test_cKO 13"
lookup_list[["E13"]][["het"]] <- "T-test_HET 13"
lookup_list[["E16"]] <- list()
lookup_list[["E16"]][["cKO"]] <- "T-test_cKO 16"
lookup_list[["E16"]][["het"]] <- "T-test_HET 16"
lookup_list[["P0"]] <- list()
lookup_list[["P0"]][["cKO"]] <- "T-test_cKO 21"
lookup_list[["P0"]][["het"]] <- "T-test_HET 21"

GO_GENE_DF <- readRDS( file = paste(base, "/other_data/GO_Term_df.RDS", sep="") )
GO_GENE_DF <- GO_GENE_DF[,c(2,1)]

eGO <- list()
dep_list <- list()
dep_count <- list()

# =======================================
# Determining the DEPs and GO enrichment 
# =======================================
for (age in c("E13", "E16", "P0")) {
  
  for (genotype in c("het", "cKO")) {
    #read data
    #Read the excel file with DEPs
    DEP_df <- read.xlsx(xlsxFile = paste(base, "other_data/proteinGroups - Reg. analysis.xlsx", sep=""), sheet = lookup_list[[age]][[genotype]], startRow = 2)
    
    dep_list[[ paste( age, genotype, "all", sep="_") ]] <- DEP_df$Genes[ which(DEP_df[,4] %in% c("down, FDR = 1%", "up, FDR = 1%")) ]
    dep_list[[ paste( age, genotype, "up", sep="_") ]] <-  DEP_df$Genes[ which(DEP_df[,4] %in% c("up, FDR = 1%")) ]
    dep_list[[ paste( age, genotype, "down", sep="_") ]] <-  DEP_df$Genes[ which(DEP_df[,4] %in% c("down, FDR = 1%")) ]

    #Separate the genes that are connected to one DEP
    gene_list = c()
    for ( gene in dep_list[[ paste( age, genotype, "all", sep="_") ]] ) {
      gene_list <- c(gene_list, strsplit(x = gene, split = ";")[[1]])
    }
    
    gene_list <- unique(gene_list)
    eGO[[ paste( age, genotype, sep="_") ]] <- enricher ( gene = gene_list, TERM2GENE = GO_GENE_DF, pvalueCutoff = 0.9, qvalueCutoff = 0.9,
                                                                   minGSSize = 10, 
                                                                   maxGSSize = 2000
                                                                )
  }
}

# write out GO analysis 
wb <- createWorkbook()
for (name in names(eGO)) {
  
  sheetName <- name
  addWorksheet(wb, sheetName)
  tmp_df <- eGO[[ sheetName ]]
  
  writeData( wb, sheetName, tmp_df )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol( tmp_df ) )
  setColWidths( wb, sheetName, cols = 1:ncol( tmp_df ), widths="auto" )
}
saveWorkbook(wb, file = paste(base, "Supplement/DEP_GO_enrichment.xlsx", sep=""), overwrite = T) 

# dep count figure
dep_list <- lapply( dep_list, function (comp){
  unique( unlist( sapply(comp, function (x) strsplit(x = x, split = ";")[[1]]) ) )
})

dep2plot <- reshape2::melt( sapply(dep_list, length ) )
dep2plot$age <- sapply( rownames(dep2plot), function (x) strsplit(x,"_")[[1]][1])
dep2plot$genotype <- sapply( rownames(dep2plot), function (x) strsplit(x,"_")[[1]][2])
dep2plot$genotype <- factor( dep2plot$genotype, levels = c("het", "cKO"))
dep2plot$dir <- sapply( rownames(dep2plot), function (x) strsplit(x,"_")[[1]][3])

ggplot(dep2plot[which(!dep2plot$dir == "all"),], aes(x=genotype, y=value, fill = dir)) + geom_bar(stat="identity") + 
  facet_wrap(~age) + theme_classic() + scale_fill_manual( values = c("up" = "goldenrod", "down" = "cornflowerblue")) + 
  ggtitle("#DEPs FDR 1%, relative to Control-MADM") + xlab("") + ylab("")
ggsave  ( filename = paste(base, "/plots/Main/Fig3C_DEP_number_bar_plot.pdf", sep=""), width=4, height=4)

# ====================================
# plot overlap cKO / HET for each age
# ====================================
plot_list <- list()
for (age in c("E13", "E16", "P0")) {
  
  tmp_deg_list <- list(c(dep_list[[ paste(age, "cKO_all", sep="_") ]]), c(dep_list[[ paste(age, "het_all", sep="_") ]]) )
  names(tmp_deg_list) <- c(paste(age, "cKO", sep="_"), paste(age, "het", sep="_") )
  plot_list[[ age ]] <- ggVennDiagram( x = tmp_deg_list, label_geom = "label" ) + 
    scale_fill_gradient(low="white",high = "white") + theme(legend.position = "NONE") + coord_cartesian()
  
}

venn_plot <- plot_grid( plotlist = plot_list, ncol = 1, greedy = F )
ggsave  ( filename = paste(base, "/plots/Supplement/FigS3F_DEP_venn_plot_comp.pdf", sep=""), plot = venn_plot, width=6, height=20)

# ===================================================
# plot overlap for ages separately for cKO / HET
# ===================================================
plot_list <- list()
for (comp in c("cKO", "het")) {
  
  tmp_deg_list <- list(c(dep_list[[ paste("E13", comp, "all", sep="_") ]]), 
                       c(dep_list[[ paste("E16", comp, "all", sep="_") ]]), 
                       c(dep_list[[ paste("P0", comp, "all", sep="_") ]]) )
  
  names(tmp_deg_list) <- c("E13", "E16", "P0")
  
  plot_list[[ comp ]] <- ggVennDiagram( x = tmp_deg_list, label_geom = "label" ) + 
    scale_fill_gradient(low="white",high = "white") + theme(legend.position = "NONE") + coord_cartesian() +
    ggtitle( comp )
  
}

venn_plot <- plot_grid( plotlist = plot_list, ncol = 1, greedy = F )
ggsave  ( filename = paste(base, "/plots/Main/Fig3E_DEP_venn_plot_ages.pdf", sep=""), plot = venn_plot, width=3, height=5)

# write out number of DEPs
wb <- createWorkbook()

  sheetName <- "DEP number"
  addWorksheet(wb, sheetName)

  writeData( wb, sheetName, dep2plot )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol( dep2plot ) )
  setColWidths( wb, sheetName, cols = 1:ncol( dep2plot ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/DEP_numbers.xlsx", sep=""), overwrite = T) 

# ==================
# GO term plots
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

# write out GO classification
wb <- createWorkbook()

sheetName <- "GO groups"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, df )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( df ) )
setColWidths( wb, sheetName, cols = 1:ncol( df ), widths="auto" )

saveWorkbook(wb, file = paste(base, "Supplement/DEP_GO_classification.xlsx", sep=""), overwrite = T) 

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
GOgroup_df2plot$GO_group <- factor( GOgroup_df2plot$GO_group, levels = rev( c( "metabol", "mitoch", "cell_cycle", "axon", "apop" ) ) )
GOgroup_df2plot$genotype <- factor( GOgroup_df2plot$genotype, levels = c("het", "cKO"))

ggplot(GOgroup_df2plot, aes(x=genotype, y=GO_group, fill=value)) + geom_tile() + 
  # scale_fill_gradient2(low="blue", mid="white", high = "gold", midpoint = 1, limits = c(0,2.5)) + 
  scale_fill_gradient(low="white", high = "black", limits = c(0,2.5)) + 
  facet_wrap(~age) + theme_classic() + ggtitle("transcriptomic GO term analysis")
ggsave  ( filename = paste(base, "/plots/Supplement/FigS3G_DEP_GO_score_plot.pdf", sep=""), width=4, height=4)


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
#   [1] ggVennDiagram_1.5.2    cowplot_1.1.3          clusterProfiler_4.10.1 ggplot2_3.5.0          openxlsx_4.2.5.2      
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               bitops_1.0-7            gson_0.1.0              shadowtext_0.1.3       
# [5] gridExtra_2.3           rlang_1.1.3             magrittr_2.0.3          DOSE_3.28.2            
# [9] compiler_4.3.2          RSQLite_2.3.6           systemfonts_1.0.6       png_0.1-8              
# [13] vctrs_0.6.5             reshape2_1.4.4          stringr_1.5.1           pkgconfig_2.0.3        
# [17] crayon_1.5.2            fastmap_1.1.1           XVector_0.42.0          labeling_0.4.3         
# [21] ggraph_2.2.1            utf8_1.2.4              HDO.db_0.99.1           enrichplot_1.22.0      
# [25] ragg_1.3.0              purrr_1.0.2             bit_4.0.5               zlibbioc_1.48.2        
# [29] cachem_1.0.8            aplot_0.2.2             jsonlite_1.8.8          GenomeInfoDb_1.38.8    
# [33] blob_1.2.4              BiocParallel_1.36.0     tweenr_2.0.3            parallel_4.3.2         
# [37] R6_2.5.1                stringi_1.8.3           RColorBrewer_1.1-3      GOSemSim_2.28.1        
# [41] Rcpp_1.0.12             IRanges_2.36.0          Matrix_1.6-5            splines_4.3.2          
# [45] igraph_2.1.4            tidyselect_1.2.1        qvalue_2.34.0           rstudioapi_0.16.0      
# [49] viridis_0.6.5           codetools_0.2-19        lattice_0.21-9          tibble_3.2.1           
# [53] plyr_1.8.9              treeio_1.26.0           Biobase_2.62.0          withr_3.0.0            
# [57] KEGGREST_1.42.0         gridGraphics_0.5-1      scatterpie_0.2.2        polyclip_1.10-6        
# [61] zip_2.3.1               Biostrings_2.70.3       pillar_1.9.0            ggtree_3.10.1          
# [65] packrat_0.9.2           stats4_4.3.2            ggfun_0.1.4             generics_0.1.3         
# [69] RCurl_1.98-1.14         S4Vectors_0.40.2        munsell_0.5.0           scales_1.3.0           
# [73] tidytree_0.4.6          glue_1.7.0              lazyeval_0.2.2          tools_4.3.2            
# [77] data.table_1.15.0       fgsea_1.28.0            fs_1.6.3                graphlayouts_1.1.1     
# [81] fastmatch_1.1-4         tidygraph_1.3.1         grid_4.3.2              tidyr_1.3.1            
# [85] ape_5.8                 AnnotationDbi_1.64.1    colorspace_2.1-0        nlme_3.1-163           
# [89] GenomeInfoDbData_1.2.11 patchwork_1.2.0         ggforce_0.4.2           cli_3.6.2              
# [93] textshaping_0.3.7       fansi_1.0.6             viridisLite_0.4.2       dplyr_1.1.4            
# [97] gtable_0.3.4            yulab.utils_0.1.4       digest_0.6.34           BiocGenerics_0.48.1    
# [101] ggrepel_0.9.5           ggplotify_0.1.2         farver_2.1.1            memoise_2.0.1          
# [105] lifecycle_1.0.4         httr_1.4.7              GO.db_3.18.0            bit64_4.0.5            
# [109] MASS_7.3-60      
