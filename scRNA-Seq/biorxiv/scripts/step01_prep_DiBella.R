# ==============================================================
# first step of the scRNA-Seq analysis
# script to prepare DiBella data for downstream analyses
# can be used for all ages or restricted to specific data sets
# ==============================================================

library (Seurat)
library (data.table)

# =======================
# set up I/O folders
# =======================
single_cell_base <- "./scRNA-Seq/biorxiv/"

# =========================================================
# read meta data 
# metaData_scDevSC.txt downloaded from single cell portal
# =========================================================
# at the time of download (03/2022):
# https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
meta_data <- fread(paste(single_cell_base, "DiBella/metaData_scDevSC.txt", sep=""))
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[2:nrow(meta_data),]

# remove spaces from ttype IDs - makes problems with saving intermediate files
meta_data$New_cellType <- gsub( pattern = " ", replacement = "_", x = meta_data$New_cellType)

# ============================================
# link gene count files to developmental age 
# ============================================

# files come from GSE153162_RAW.tar downloaded from GEO
file_list <- c("GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5", "GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5",
               "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5", "GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5",
               "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5", "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
                "GSM5277844_E17_5_filtered_feature_bc_matrix.h5", "GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5",
               "GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5", "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
               "GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5", "GSM5277843_E10_v1_filtered_feature_bc_matrix.h5" )

age_list <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18_S1", "E18_S3", "P1_S1", "P1", "E10")
ct_prefix <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E18", "P1", "P1", "E10")

seurat.obj.lst <- list()

# ===================================
# create Seurat object from raw data
# ===================================

# focus on P1 here
for ( i in c( 10, 11 )) {  

  cell_types_to_remove <- c( "Low_quality_cells", "Doublet" )
  
  fileName <- file_list[i]
  
  sc.data <- Read10X_h5(filename = paste(base, "DiBella/", fileName, sep=""))

  sub_meta_data <- meta_data[ which(meta_data$biosample_id ==  age_list[i] ), ]
	
  # remove cell types with less than 50 cells
  cell_types_to_remove <- c( cell_types_to_remove, names ( which( table(sub_meta_data$New_cellType) < 50 ) ) )
  sub_meta_data <- sub_meta_data[ which(!sub_meta_data$New_cellType %in% cell_types_to_remove), ]

  cellIDs <- sub_meta_data$NAME
	
  cellIDs <- gsub( pattern = "-1", replacement = "", x = cellIDs )
	
  cellIDs <- sapply(cellIDs, function (x) {
    tmp <- strsplit(x = x, split = "_")[[1]]
    return ( tmp[length(tmp)] )
  })
  
  cellIDs <- paste( cellIDs, 1, sep="-" )
  
  rownames(sub_meta_data) <- cellIDs
	
  colIDX <- which( colnames(sc.data) %in% cellIDs)
  sc.data <- sc.data[,colIDX]
  sub_meta_data <- sub_meta_data[colnames(sc.data),]
	
  sub_meta_data$New_cellType <- paste( sub_meta_data$New_cellType, ct_prefix[i], sep="_")
  sub_meta_data$age <- ct_prefix[i]
    
  message ( paste( age_list[i], fileName, nrow( sub_meta_data ), ncol(sc.data), nrow(sc.data), ct_prefix[i] ) )
	
  seurat.obj.lst[[ age_list[i] ]] <-  CreateSeuratObject(counts = sc.data, project = age_list[i], min.cells = 3, min.features = 200)
  seurat.obj.lst[[ age_list[i] ]] <- AddMetaData( object = seurat.obj.lst[[ age_list[i] ]], metadata = sub_meta_data[, c("age", "New_cellType", "Phase")] )
	
}

# not provided on GitHub
saveRDS( object = seurat.obj.lst, file = paste(base, "RDS_files/DiBella.seurat.obj.list.RDS", sep=""))

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
#   [1] data.table_1.15.0  Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-3          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.3            magrittr_2.0.3        
# [6] RcppAnnoy_0.0.22       spatstat.geom_3.2-8    matrixStats_1.2.0      ggridges_0.5.6         compiler_4.3.2        
# [11] png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         hdf5r_1.3.10           stringr_1.5.1         
# [16] pkgconfig_2.0.3        fastmap_1.1.1          ellipsis_0.3.2         utf8_1.2.4             promises_1.2.1        
# [21] bit_4.0.5              purrr_1.0.2            jsonlite_1.8.8         goftest_1.2-3          later_1.3.2           
# [26] spatstat.utils_3.1-0   irlba_2.3.5.1          parallel_4.3.2         cluster_2.1.4          R6_2.5.1              
# [31] ica_1.0-3              stringi_1.8.3          RColorBrewer_1.1-3     spatstat.data_3.0-4    reticulate_1.35.0     
# [36] parallelly_1.37.0      lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.12            tensor_1.5            
# [41] future.apply_1.11.1    zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.14          Matrix_1.6-5          
# [46] splines_4.3.2          igraph_2.0.2           tidyselect_1.2.1       rstudioapi_0.16.0      abind_1.4-5           
# [51] spatstat.random_3.2-2  codetools_0.2-19       miniUI_0.1.1.1         spatstat.explore_3.2-6 listenv_0.9.1         
# [56] lattice_0.21-9         tibble_3.2.1           plyr_1.8.9             shiny_1.8.0            ROCR_1.0-11           
# [61] Rtsne_0.17             future_1.33.1          fastDummies_1.7.3      survival_3.5-7         polyclip_1.10-6       
# [66] fitdistrplus_1.1-11    pillar_1.9.0           packrat_0.9.2          KernSmooth_2.23-22     plotly_4.10.4         
# [71] generics_0.1.3         RcppHNSW_0.6.0         ggplot2_3.5.0          munsell_0.5.0          scales_1.3.0          
# [76] globals_0.16.2         xtable_1.8-4           glue_1.7.0             lazyeval_0.2.2         tools_4.3.2           
# [81] RSpectra_0.16-1        RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3         
# [86] grid_4.3.2             tidyr_1.3.1            colorspace_2.1-0       nlme_3.1-163           patchwork_1.2.0       
# [91] cli_3.6.2              spatstat.sparse_3.0-3  spam_2.10-0            fansi_1.0.6            viridisLite_0.4.2     
# [96] dplyr_1.1.4            uwot_0.1.16            gtable_0.3.4           digest_0.6.34          progressr_0.14.0      
# [101] ggrepel_0.9.5          htmlwidgets_1.6.4      htmltools_0.5.7        lifecycle_1.0.4        httr_1.4.7            
# [106] mime_0.12              bit64_4.0.5            MASS_7.3-60    
