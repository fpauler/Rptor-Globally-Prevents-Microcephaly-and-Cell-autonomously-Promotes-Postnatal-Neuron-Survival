# ===================================================================
# second step of the analysis
# identifies cell type marker genes from DiBella et al. Nature 2021
# will be used for ROS analysis
# ===================================================================

library(Seurat)
library (ggplot2)

# =======================
# set up I/O folders
# =======================
base <- "./bulk_RNASeq/biorxiv/"
scRNA_base <- "./scRNA-Seq/biorxiv/"

# ===================
# prepare reference
# ===================
# Note: this file comes from the single cell analysis: step01_prep_DiBella.R
#       this file is >50Mb and is not provided on GitHub
seurat.obj.lst <- readRDS( file = paste( scRNA_base, "RDS_files/DiBella.seurat.obj.list.RDS", sep="") )

# NOTE: I only use 1 P1 data set as reference
reference <- seurat.obj.lst[[1]]

# standard workflow and some QC
reference <- NormalizeData(reference, normalization.method = "LogNormalize", scale.factor = 10000)
reference <- FindVariableFeatures( reference, selection.method = "vst", nfeatures = 2000 )
reference <- ScaleData( reference )
reference <- RunPCA( reference )
reference <- RunUMAP( reference, dims = 1:15 )

DimPlot( reference, reduction = "umap", group.by = "orig.ident")
ggsave( paste(base, "plots/QC/Reference_P1_umap.pdf", sep=""))
DimPlot( reference, reduction = "umap", group.by = "New_cellType", label=T)
ggsave( paste(base, "plots/QC/Reference_P1_umap_New_cellType.pdf", sep=""), width=6, height=5)

# DL cells separate into multiple types, combine here
reference@meta.data$mod_celltype <- reference@meta.data$New_cellType
reference@meta.data$mod_celltype[ which( reference@meta.data$mod_celltype %in% c("CThPN_P1", "SCPN_P1", "Layer_4_P1", "Layer_6b_P1", "NP_P1" ) )] <- "DL_CPN_P1"

DimPlot( reference, reduction = "umap", group.by = "mod_celltype", label=T)
ggsave( paste(base, "plots/QC/Reference_P1_umap_mod_celltype.pdf", sep=""), width=6, height=5)

Idents(reference) <- "mod_celltype"
markers <- FindAllMarkers( object = reference, only.pos = T)

marker_list <- list()
marker_df <- data.frame()
for (cell_type in unique(markers$cluster)) {
  marker_list[[ cell_type ]] <- markers$gene[which(markers$cluster == cell_type )[1:200]]
  marker_df <- rbind( marker_df, data.frame( term = cell_type, gene = marker_list[[ cell_type ]] ))
}

# check for specificity of markers
table( table( marker_df$gene) )
#    1    2    3    4 
# 1113  137    3    1 

saveRDS( object = marker_df, file = paste( base, "other_data/DiBella.P1.marker_df.RDS", sep=""))

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
#   [1] ggplot2_3.5.0      Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.3            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.22       spatstat.geom_3.2-8    matrixStats_1.2.0      ggridges_0.5.6         compiler_4.3.2        
# [11] systemfonts_1.0.6      png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1         
# [16] pkgconfig_2.0.3        fastmap_1.1.1          ellipsis_0.3.2         labeling_0.4.3         utf8_1.2.4            
# [21] promises_1.2.1         ragg_1.3.0             purrr_1.0.2            jsonlite_1.8.8         goftest_1.2-3         
# [26] later_1.3.2            spatstat.utils_3.1-0   irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.4         
# [31] R6_2.5.1               ica_1.0-3              stringi_1.8.3          RColorBrewer_1.1-3     spatstat.data_3.0-4   
# [36] limma_3.58.1           reticulate_1.35.0      parallelly_1.37.0      lmtest_0.9-40          scattermore_1.2       
# [41] Rcpp_1.0.12            tensor_1.5             future.apply_1.11.1    zoo_1.8-12             sctransform_0.4.1     
# [46] httpuv_1.6.14          Matrix_1.6-5           splines_4.3.2          igraph_2.0.2           tidyselect_1.2.1      
# [51] rstudioapi_0.16.0      abind_1.4-5            spatstat.random_3.2-2  codetools_0.2-19       miniUI_0.1.1.1        
# [56] spatstat.explore_3.2-6 listenv_0.9.1          lattice_0.21-9         tibble_3.2.1           plyr_1.8.9            
# [61] withr_3.0.0            shiny_1.8.0            ROCR_1.0-11            Rtsne_0.17             future_1.33.1         
# [66] fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6        fitdistrplus_1.1-11    pillar_1.9.0          
# [71] packrat_0.9.2          KernSmooth_2.23-22     plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
# [76] munsell_0.5.0          scales_1.3.0           globals_0.16.2         xtable_1.8-4           glue_1.7.0            
# [81] lazyeval_0.2.2         tools_4.3.2            data.table_1.15.0      RSpectra_0.16-1        RANN_2.6.1            
# [86] leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3          grid_4.3.2             tidyr_1.3.1           
# [91] colorspace_2.1-0       nlme_3.1-163           patchwork_1.2.0        presto_1.0.0           cli_3.6.2             
# [96] spatstat.sparse_3.0-3  textshaping_0.3.7      spam_2.10-0            fansi_1.0.6            viridisLite_0.4.2     
# [101] dplyr_1.1.4            uwot_0.1.16            gtable_0.3.4           digest_0.6.34          progressr_0.14.0      
# [106] ggrepel_0.9.5          farver_2.1.1           htmlwidgets_1.6.4      htmltools_0.5.7        lifecycle_1.0.4       
# [111] httr_1.4.7             statmod_1.5.0          mime_0.12              MASS_7.3-60 
