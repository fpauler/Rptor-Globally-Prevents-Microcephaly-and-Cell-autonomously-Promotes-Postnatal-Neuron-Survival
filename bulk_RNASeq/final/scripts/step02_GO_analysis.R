# ========================================================
# step 2: of revision re-analysis
# contrast DEGs and DEPs GO groups
# prepare for STRING analysis
# ========================================================

library (openxlsx)
library (clusterProfiler)
library (ggplot2)
library (tidyr)
library (dplyr)

# =======================
# set up I/O folders
# =======================
bulk_base <- "./bulk_RNASeq/biorxiv/" # initial analysis
main_base <- "./bulk_RNASeq/final/" # main working folder

# ==========================
# read GO terms from step1
# ==========================
GO_GENE_DF <- readRDS(  paste(main_base, "/RDS_files/GO_Term_df_collapsed.RDS", sep="")  )
colnames( GO_GENE_DF ) <- c("TERM", "GENE")
# give some group a simpler name
GO_GENE_DF$TERM[ which(GO_GENE_DF$TERM == "reactive oxygen")] <- "ROS"
GO_GENE_DF$TERM[ which(GO_GENE_DF$TERM == "epidermal growth factor")] <- "Egfr"
GO_GENE_DF$TERM[ which(GO_GENE_DF$TERM == "mitochondr")] <- "mito"

# ===========================
# container for GO analysis
# ===========================
eGO <- list()

# ========================
# Transcriptomic analysis
# ========================

deg_list <- list()
for (age in c("E13", "E16", "P0")) {
  for ( group in c("cKO", "HET") ) {
    deg_df <- read.xlsx( paste(bulk_base, "Supplement/DEG_raw_DEG_stats.xlsx", sep=""), sheet = paste( age, group, "WT", sep="_" ))
    
    deg_list[[paste(age, toupper(group), sep="_")]] <- deg_df$symbol[ which(deg_df$padj < 0.05) ]

  }
}

# combine DEGs from different developmental stages, separately for cKO and Het
gene_list_cKO <- unique( unlist( deg_list[ c("E13_CKO", "E16_CKO", "P0_CKO") ] ) )
gene_list_het <- unique( unlist( deg_list[ c("E13_HET", "E16_HET", "P0_HET") ] ) )

dgl <- list()
dgl[["common"]] <- intersect( gene_list_cKO, gene_list_het)
dgl[["het"]] <- setdiff( gene_list_het, gene_list_cKO)
dgl[["cKO"]] <- setdiff( gene_list_cKO, gene_list_het )

for (comp in c("common", "het", "cKO")) {
  id <- paste("deg", comp, sep="_")
  tmp_eGO <- enricher ( gene = dgl[[ comp ]], 
                        TERM2GENE = GO_GENE_DF[,c("TERM", "GENE")], 
                        pvalueCutoff = 0.99, qvalueCutoff = 0.99, 
                        pAdjustMethod = "BH",
                        minGSSize = 10, 
                        maxGSSize = 3000 )
  
  eGO[[ id ]] <- tmp_eGO@result
  
}

# ===================
# Protein analysis
# ===================

# define lookup list to where to find the dataset for each age and genotype
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

dep_list <- list()
for (age in c("E13", "E16", "P0")) {
  
  for (genotype in c("het", "cKO")) {
    # read data
    # Read the excel file with DEPs
    DEP_df <- read.xlsx(xlsxFile = paste(bulk_base, "/other_data/proteinGroups - Reg. analysis.xlsx", sep=""), sheet = lookup_list[[age]][[genotype]], startRow = 2)
    
    dep_vec <- DEP_df$Genes[ which(DEP_df[,4] %in% c("down, FDR = 1%", "up, FDR = 1%")) ]

    # Separate the genes that are connected to one DEP
    gene_list = c()
    for ( gene in dep_vec ) {
      gene_list <- c(gene_list, strsplit(x = gene, split = ";")[[1]])
    }
    
    dep_list[[paste(age, toupper(genotype), sep="_")]] <- unique( gene_list )
    
  }
}

dep_list_cKO <- unique( unlist( dep_list[ c("E13_CKO", "E16_CKO", "P0_CKO") ] ) )
dep_list_het <- unique( unlist( dep_list[ c("E13_HET", "E16_HET", "P0_HET") ] ) )

dep_list <- list()
dep_list[["common"]] <- intersect( dep_list_cKO, dep_list_het)
dep_list[["het"]] <- setdiff( dep_list_het, dep_list_cKO)
dep_list[["cKO"]] <- setdiff( dep_list_cKO, dep_list_het )

for (comp in c("common", "het", "cKO")) {
  id <- paste("dep", comp, sep="_")
  tmp_eGO <- enricher ( gene = dep_list[[ comp ]], 
                      TERM2GENE = GO_GENE_DF[,c("TERM", "GENE")], 
                      pvalueCutoff = 0.99, qvalueCutoff = 0.99, 
                      pAdjustMethod = "BH",
                      minGSSize = 10, 
                      maxGSSize = 3000 )
  eGO[[ id ]] <- tmp_eGO@result
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
saveWorkbook(wb, file = paste(main_base, "/Supplement/GO_enrichment.xlsx", sep=""), overwrite = T) 

# plot fractions of genes
df2plot <- data.frame()
for (name in names(eGO)) {
  
  tmp_df <- eGO[[ name ]]
  
  type <- strsplit( name, "_")[[1]][1]
  genotype <- strsplit( name, "_")[[1]][2]
  
  fractions <- as.numeric( sapply( tmp_df$GeneRatio, function(z) eval(parse(text = z))) )
  
  df2plot <- rbind( df2plot, data.frame( Description=tmp_df$Description, type = type, genotype = genotype, frac = fractions ) )
  
}

# add 0 for missing values
for (descr in unique( df2plot$Description) ) {
  for (type in unique( df2plot$type)) {
    for (gen in unique( df2plot$genotype)) {
      idx <- which( df2plot$Description == descr & df2plot$type == type & df2plot$genotype == gen )
      if (length(idx) == 0) {
        message( descr, type, gen)
        df2plot <- rbind( df2plot, data.frame(Description = descr, type = type, genotype = gen, frac = 0))
      }
    }
  }
}

sig_pw <- c("Egfr", "Notch", "TOR/PI3K", "JAK-STAT", "MAPK", "p53")

df2plot_sig <- df2plot[which(df2plot$Description %in% sig_pw),]

sum_frac <- df2plot_sig %>% group_by( c(Description)) %>% summarise( sum(frac))
colnames(sum_frac) <- c("Description", "sum")
df2plot_sig$Description <- factor( df2plot_sig$Description, levels = sum_frac$Description[ order(sum_frac$sum, decreasing = F) ] )

ggplot( df2plot_sig, aes(x=frac, y=Description, fill=genotype)) + geom_bar( stat="identity", position="dodge") +
  facet_grid( ~ type) + theme_classic() + ggtitle("fraction of genes in category")
ggsave( filename = paste(main_base, "/plots/Main/Fig3e_GO_gene_frac_sig.pdf", sep=""), width=6, height=5 )

df2plot_BP <- df2plot[ which( !df2plot$Description %in% sig_pw ) , ]

sum_frac <- df2plot_BP %>% group_by( c(Description)) %>% summarise( sum(frac))
colnames(sum_frac) <- c("Description", "sum")
df2plot_BP$Description <- factor( df2plot_BP$Description, levels = sum_frac$Description[ order(sum_frac$sum, decreasing = F) ] )

ggplot( df2plot_BP, aes(x=frac, y=Description, fill=genotype)) + geom_bar( stat="identity", position="dodge") +
  facet_grid( ~ type) + theme_classic() + ggtitle("fraction of genes in category")
ggsave( filename = paste(main_base, "/plots/Main/Fig3d_GO_gene_frac_BP.pdf", sep=""), width=6, height=5 )

# =============================================================================================
# combine DEGs and DEPs for an holistic view of the changes in Rptor-MADM and cKO-Rptor-MADM
# note that multiple GO group associations are combied to a comma separated list
# this list is used in a later step for STRING network analysis
# =============================================================================================

combined_genes <- eGO[c("deg_common", "deg_het", "dep_common", "dep_het")] %>%
  lapply(function(df) select(df, ID, geneID)) %>%
  bind_rows() %>%
  separate_rows(geneID, sep = "/") %>%
  rename(gene = geneID, Description = ID) %>%
  distinct(gene, Description) %>%
  group_by(gene) %>%
  summarise(
    Description = paste(unique(Description), collapse = ", "),
    .groups = "drop"
  )

write.csv( x = combined_genes, file = paste(main_base, "/Supplement/het_genes_4nw.csv", sep="") )

combined_genes <- eGO[c("deg_common", "deg_cKO", "dep_common", "dep_cKO")] %>%
  lapply(function(df) select(df, ID, geneID)) %>%
  bind_rows() %>%
  separate_rows(geneID, sep = "/") %>%
  rename(gene = geneID, Description = ID) %>%
  distinct(gene, Description) %>%
  group_by(gene) %>%
  summarise(
    Description = paste(unique(Description), collapse = ", "),
    .groups = "drop"
  )

write.csv( x = combined_genes, file = paste(main_base, "/Supplement/cKO_genes_4nw.csv", sep="") )

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
#   [1] dplyr_1.1.4            tidyr_1.3.1            ggplot2_3.5.2          clusterProfiler_4.10.1 openxlsx_4.2.5.2      
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               bitops_1.0-7            gson_0.1.0              shadowtext_0.1.3        gridExtra_2.3          
# [6] rlang_1.1.3             magrittr_2.0.3          DOSE_3.28.2             compiler_4.3.2          RSQLite_2.3.6          
# [11] systemfonts_1.0.6       png_0.1-8               vctrs_0.6.5             reshape2_1.4.4          stringr_1.5.1          
# [16] pkgconfig_2.0.3         crayon_1.5.2            fastmap_1.1.1           XVector_0.42.0          labeling_0.4.3         
# [21] ggraph_2.2.1            utf8_1.2.4              HDO.db_0.99.1           enrichplot_1.22.0       ragg_1.3.0             
# [26] purrr_1.0.2             bit_4.0.5               zlibbioc_1.48.2         cachem_1.0.8            aplot_0.2.2            
# [31] jsonlite_1.8.8          GenomeInfoDb_1.38.8     blob_1.2.4              BiocParallel_1.36.0     tweenr_2.0.3           
# [36] parallel_4.3.2          R6_2.5.1                stringi_1.8.3           RColorBrewer_1.1-3      pkgload_1.3.4          
# [41] GOSemSim_2.28.1         Rcpp_1.0.12             IRanges_2.36.0          Matrix_1.6-5            splines_4.3.2          
# [46] igraph_2.1.4            tidyselect_1.2.1        qvalue_2.34.0           rstudioapi_0.16.0       viridis_0.6.5          
# [51] codetools_0.2-19        lattice_0.21-9          tibble_3.2.1            plyr_1.8.9              Biobase_2.62.0         
# [56] treeio_1.26.0           withr_3.0.0             KEGGREST_1.42.0         gridGraphics_0.5-1      scatterpie_0.2.2       
# [61] polyclip_1.10-6         zip_2.3.1               Biostrings_2.70.3       pillar_1.9.0            ggtree_3.10.1          
# [66] packrat_0.9.2           stats4_4.3.2            ggfun_0.1.4             generics_0.1.3          RCurl_1.98-1.14        
# [71] S4Vectors_0.40.2        scales_1.4.0            tidytree_0.4.6          glue_1.7.0              lazyeval_0.2.2         
# [76] tools_4.3.2             data.table_1.15.0       fgsea_1.28.0            fs_1.6.3                graphlayouts_1.1.1     
# [81] fastmatch_1.1-4         tidygraph_1.3.1         cowplot_1.1.3           grid_4.3.2              ape_5.8                
# [86] AnnotationDbi_1.64.1    colorspace_2.1-0        nlme_3.1-163            GenomeInfoDbData_1.2.11 patchwork_1.2.0        
# [91] ggforce_0.4.2           cli_3.6.2               textshaping_0.3.7       fansi_1.0.6             viridisLite_0.4.2      
# [96] gtable_0.3.6            yulab.utils_0.1.4       digest_0.6.34           BiocGenerics_0.48.1     ggrepel_0.9.5          
# [101] ggplotify_0.1.2         farver_2.1.1            memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7             
# [106] GO.db_3.18.0            bit64_4.0.5             MASS_7.3-60  
