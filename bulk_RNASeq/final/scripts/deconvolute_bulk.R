# ========================================================================================
# script to deconvolute bulk-seq data using diBella et al Nature 2021 as a reference
# Note: this analysis was done using R version 4.5.3
# ========================================================================================

library (Seurat)
library (data.table)
library (ggplot2)
library (openxlsx)
library (BayesPrism)
library (zCompositions)
library (compositions)
library (limma)
library (tidyr)
library (dplyr)
library (tibble)
library (DESeq2)

# ==========
# functions 
# ==========

## testing for cell type proportion changes
## ============================================================================
## CoDA-based differential abundance testing for deconvolution-derived
## cell-type proportions
##
## Input:
##   prop_mat  - samples x cell-types numeric matrix/data.frame, rows should
##               sum to (approximately) 1. Rownames = sample IDs,
##               colnames = cell-type names.
##   group     - factor/character vector of group labels, length == nrow(prop_mat),
##               in the same sample order as prop_mat.
##
## Output:
##   - MANOVA result on ilr coordinates (global compositional shift test)
##   - Per-cell-type results table (CLR + limma/eBayes moderated test)
##   - A boxplot of clr-transformed proportions per cell type, faceted by group
##   - A PCA (compositional biplot) of samples colored by group
## ============================================================================

## ---- 1. Main function ---------------------------------------------------

#' Run CoDA differential abundance testing on cell-type proportions
#'
#' @param prop_mat  samples x cell-types matrix/data.frame of proportions
#' @param group     vector of group labels, same length/order as rows of prop_mat
#' @param zero_method  "CZM" (default) passed to zCompositions::cmultRepl, or
#'                      "pseudocount" for a simple additive pseudocount fallback
#' @param p_adjust_method  multiple testing correction for per-cell-type table
#'
#' @return a list with elements
#' if more than 2 cell types:
#'   manova        - summary.manova object (Pillai's trace test) on ilr coords
#'   results_table  - data.frame: cell_type, logFC, AveExpr, t, P.Value, adj.P.Val
#'   clr_coords     - matrix of clr-transformed proportions (for custom plotting)
#'   ilr_coords     - matrix of ilr-transformed proportions
#'   boxplot        - ggplot object, clr proportions by cell type & group
#'   pca_plot       - ggplot object, compositional PCA (clr-based) colored by group
#' if 2 cell types:
#'   test        - in case of 2 group values: wilcox.test, >2 group values: kruskal.test of log ratios
#'   boxplot     - ggplot object, log ratios by cell type & group,
#'   log_values  - log ratios
#' 

run_coda_diffabundance <- function(prop_mat,
                                   group,
                                   zero_method = c("CZM", "pseudocount"),
                                   p_adjust_method = "BH", age = age) {
  
  zero_method <- match.arg(zero_method)
  prop_mat <- as.matrix(prop_mat)
  group <- factor(group)
  
  stopifnot(
    "prop_mat rownames required" = !is.null(rownames(prop_mat)),
    "prop_mat colnames required" = !is.null(colnames(prop_mat)),
    "group length must match nrow(prop_mat)" = length(group) == nrow(prop_mat),
    "need at least 2 groups" = nlevels(group) >= 2
  )
  
  # Given all this, for exactly 2 cell types, the ilr/clr machinery is mostly ceremony — 
  # I just compute a single log-ratio and test that:
  if (ncol(prop_mat) == 2){
    
    log_ratio <- log(prop_mat[, 1] / prop_mat[, 2])
    
    if (nlevels(group) == 2) {
      test <- wilcox.test(log_ratio ~ group)
    } else {
      test <- kruskal.test(log_ratio ~ group)
    }
    
    log_df <- data.frame( log_ratio = log_ratio, group = group)
    
    boxplot_gg <- ggplot(log_df, aes(x = group, y = log_ratio, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
      labs(
        title = paste(age, paste(colnames(prop_mat), collapse=" vs "), sep = " "),
        x = NULL, y = "log-transformed proportion"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    list(
      test        = test,
      boxplot     = boxplot_gg,
      log_values  = log_df
    ) -> return_list
    
  } else {
    
    # Re-close rows to sum to 1 in case of rounding drift
    prop_mat <- prop_mat / rowSums(prop_mat)
    
    ## ---- 1a. Zero replacement ----
    if (any(prop_mat == 0)) {
      if (zero_method == "CZM") {
        prop_no_zero <- cmultRepl(prop_mat, method = "CZM", label = 0, output = "prop")
      } else {
        min_nonzero <- min(prop_mat[prop_mat > 0])
        prop_no_zero <- prop_mat
        prop_no_zero[prop_no_zero == 0] <- min_nonzero / 2
        prop_no_zero <- prop_no_zero / rowSums(prop_no_zero)
      }
    } else {
      prop_no_zero <- prop_mat
    }
    
    ## ---- 1b. Transforms ----
    comp_obj <- acomp(prop_no_zero)
    
    ilr_raw <- ilr(comp_obj)
    # Force plain numeric matrix - ilr() returns an 'rmult' object (compositions
    # pkg) which can behave unpredictably inside manova()/model.matrix()
    ilr_coords <- matrix(as.numeric(ilr_raw), nrow = nrow(prop_mat),
                         ncol = ncol(ilr_raw))
    colnames(ilr_coords) <- paste0("ilr", seq_len(ncol(ilr_coords)))
    rownames(ilr_coords) <- rownames(prop_mat)
    
    clr_coords <- clr(comp_obj)
    clr_coords <- matrix(as.numeric(clr_coords), nrow = nrow(clr_coords),
                         dimnames = dimnames(clr_coords))
    rownames(clr_coords) <- rownames(prop_mat)
    
    ## ---- 1c. Global test: MANOVA on ilr coordinates ----
    manova_fit <- manova(ilr_coords ~ group)
    manova_summary <- summary(manova_fit, test = "Pillai")
    
    ## ---- 1d. Per-cell-type test: CLR + limma/eBayes ----
    design <- model.matrix(~ group)
    fit <- lmFit(t(clr_coords), design)
    fit <- eBayes(fit)
    
    # coef = last column = last group level vs reference; for >2 groups this
    # tests the first non-reference level by default. For a full multi-group
    # F-test across all group coefficients, see note below.
    if (nlevels(group) == 2) {
      results_table <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
    } else {
      # ANOVA-like F-test across all non-intercept group coefficients
      results_table <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = "F")
    }
    results_table$cell_type <- rownames(results_table)
    rownames(results_table) <- NULL
    results_table <- results_table %>% relocate(cell_type)
    
    ## ---- 1e. Boxplot: clr proportions by cell type, faceted/colored by group ----
    clr_df <- as.data.frame(clr_coords)
    clr_df$sample <- rownames(clr_df)
    clr_df$group <- group
    
    clr_long <- clr_df %>%
      pivot_longer(cols = -c(sample, group), names_to = "cell_type", values_to = "clr_value")
    
    boxplot_gg <- ggplot(clr_long, aes(x = group, y = clr_value, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
      facet_wrap(~ cell_type, scales = "free_y") +
      labs(
        title = age,
        x = NULL, y = "CLR-transformed proportion"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ## ---- 1f. Compositional PCA (biplot-style) on clr coordinates ----
    pca_res <- prcomp(clr_coords, center = TRUE, scale. = FALSE)
    var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
    
    pca_df <- as.data.frame(pca_res$x[, 1:2])
    colnames(pca_df) <- c("PC1", "PC2")
    pca_df$sample <- rownames(clr_coords)
    pca_df$group <- group
    
    pca_plot_gg <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
      geom_point(size = 3, alpha = 0.85) +
      stat_ellipse(type = "t", level = 0.68, linetype = "dashed") +
      labs(
        title = "Compositional PCA (CLR) of cell-type proportions",
        x = paste0("PC1 (", var_explained[1], "%)"),
        y = paste0("PC2 (", var_explained[2], "%)")
      ) +
      theme_bw(base_size = 11)
    
    list(
      manova        = manova_summary,
      results_table = results_table,
      clr_coords    = clr_coords,
      ilr_coords    = ilr_coords,
      boxplot       = boxplot_gg,
      pca_plot      = pca_plot_gg
    ) -> return_list
  }
  
  return( return_list )
}

# =======================
# set up I/O folders
# =======================
single_cell_base <- "./scRNA-Seq/biorxiv/"
bulk_base <- "./bulk_RNASeq/biorxiv/" # initial analysis
main_base <- "./bulk_RNASeq/final/" # main working folder

# read a set of genes that I found to work with BayesPrism
gene_list <- readRDS( file = paste( main_base, "other_data/gene_list.RDS", sep=" ))

# =============================
# prepare the bulk count data
# =============================

#prepare environment hash table as lookup for ENSMUSG - Symbol
#this comes from STAR index folder and was prepared solely by using info from Gencode
#initialise variables
ensmusg_symbol <- new.env()
ensmusg_chr <- new.env()
ensmusg_entrezid <- new.env()

# read conversion table - this part of the analysis was performed with:
# genome sequence (GRCm38.p5), gene annotation conversion tables and Gencode M16

ensmusg_symbol_chr <- read.table(paste(bulk_base, "/annotation/exchange.table_M16.tsv", sep=""), fill = T, header = F, stringsAsFactors = F)
colnames(ensmusg_symbol_chr) <- c("ENSMUSG", "ENSMUST", "SYMBOL", "chr", "ENTREZID")

#fill conversion table - handled just like a list but apparently much faster
for (x in 1:nrow(ensmusg_symbol_chr)) {
  #deal here with multi mappers - one SYMBOL mapped to many ENSMUSG
  #a comma in the SYMBOL name indicates these multi-mappers
  if (is.null(ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]])) {
    ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] <- ensmusg_symbol_chr$SYMBOL[x]
  } else {
    
    if (!ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] == ensmusg_symbol_chr$SYMBOL[x]) {
      ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]] <- paste(ensmusg_symbol[[ensmusg_symbol_chr$ENSMUSG[x]]], ensmusg_symbol_chr$SYMBOL[x], sep=",")
    }
  }
}

ensmusg_symbol[["Rptor_del"]] <- "Rptor_del"

# ==========================
# read the gene count data
# ==========================

# Folder with gene count files - give full path
gene_count_base=paste(bulk_base, "/Rptor_SMARTer.files/", sep="")

# EXEL file with meta data - give full path
fn_sample_list=paste(bulk_base, "/Supplement/Rptor_sample_list.xlsx", sep="")

#read EXEL file using openxlsx
bulk_meta <- read.xlsx(xlsxFile = fn_sample_list, sheet = 1, skipEmptyRows = FALSE)
bulk_meta <- bulk_meta[ which(bulk_meta$used_for_analysis == "Y"),]

#match samples with read count files - done via sampleID
#recursive search for all files in the count folder - returns relative path
files <- list.files(path = gene_count_base, pattern="*ReadsPerGene.out.tab", recursive = T)
#make path absolute
full_fn_path_list <- lapply(files, function (x) {paste(gene_count_base, x, sep="")})
#prepare dataframe from file list - matching ID with filename
files.DF <- data.frame(fn = unlist(full_fn_path_list), stringsAsFactors = F)
files.DF$sample_id <- apply(files.DF, 1, function (x) {
  #extract sample ID from filename
  parts <- strsplit(x[1], "\\.")
  return(parts[[1]][length(parts[[1]])-3])
})

bulk_meta <- merge(bulk_meta, files.DF, by="sample_id")

#initialise dataframe holding all read counts df
all <- data.frame(stringsAsFactors = F)

#read in all counts and merge in one dataframe - by gene_id
for (line in 1:nrow(bulk_meta)) {
  this_count <- read.table(file = bulk_meta$fn[line], header = F, sep="\t", stringsAsFactors = F)
  #remove first 4 lines - only basic stats of features
  this_count<-this_count[5:nrow(this_count),c(1,2)]
  
  colnames(this_count) <- c("gene_id", line)
  
  if (nrow(all) == 0) {
    all <- this_count
  } else {
    all<-merge(all, this_count, by="gene_id")
  }
}

#give meaningful column names - sample_id in this case
colnames(all)<-c("gene_id", bulk_meta$sample_id)

# exchange ensembl IDs with symbols here
# note that duplicated symbol entries are removed!
gn <- all$gene_id
symbol <- sapply(gn, function (x) ensmusg_symbol[[x]])
symbol <- as.character( symbol )
gn <- data.frame( ensembl = gn, symbol = symbol)
gn <- gn[!(duplicated(gn$symbol) | duplicated(gn$symbol, fromLast = TRUE)), ]
rownames(gn) <- gn$ensembl

rownames(all) <- all$gene_id

all<-all[gn$ensembl,]
rownames(all) <- gn$symbol

# =================================
# prepare the scRNA-Seq reference
# =================================

# metaData_scDevSC.txt downloaded from single cell portal
# at the time of download (03/2022):
# https://singlecell.broadinstitute.org/single_cell/study/SCP1290/molecular-logic-of-cellular-diversification-in-the-mammalian-cerebral-cortex
meta_data <- fread( paste( single_cell_base, "DiBella/metaData_scDevSC.txt", sep="") )
meta_data <- as.data.frame(meta_data)
meta_data <- meta_data[2:nrow(meta_data),]

# remove spaces from ttype IDs - makes problems with saving intermediate files
meta_data$New_cellType <- gsub( pattern = " ", replacement = "_", x = meta_data$New_cellType)

# files come from GSE153162_RAW.tar downloaded from GEO
# link al files todevelopmental age
file_list <- c("GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5", "GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5",
               "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5", "GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5",
               "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5", "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
               "GSM5277844_E17_5_filtered_feature_bc_matrix.h5", "GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5",
               "GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5", "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
               "GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5", "GSM5277843_E10_v1_filtered_feature_bc_matrix.h5" )

age_list <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18_S1", "E18_S3", "P1_S1", "P1", "E10")
ct_prefix <- c("E11", "E12", "E13", "E14", "E15", "E16", "E17", "E18", "E18", "P1", "P1", "E10")

cell_types_to_include <- c("Apical_progenitors", "Intermediate_progenitors", "Immature_neurons", 
                           "Migrating_neurons", "SCPN", "CThPN", "DL_CPN", "UL_CPN", "Layer_4", "NP", "Layer_6b",
                           "Cycling_glial_cells", "Astrocytes", "Oligodendrocytes")

# combine cell types that might be too similar to discern 
meta_data$broad_New_cellType <- meta_data$New_cellType
meta_data$broad_New_cellType[ which( meta_data$New_cellType %in% c("Apical_progenitors", "Intermediate_progenitors"))] <- "RGP/IP"
meta_data$broad_New_cellType[ which( meta_data$New_cellType %in% c("Immature_neurons", "Migrating_neurons"))] <- "Immature_neurons"
meta_data$broad_New_cellType[ which( meta_data$New_cellType %in% c("SCPN", "CThPN", "DL_CPN", "DL_CPN_1", "DL_CPN_2", "UL_CPN", "Layer_4", "NP", "Layer_6b"))] <- "Neurons"
meta_data$broad_New_cellType[ which( meta_data$New_cellType %in% c("Astrocytes", "Oligodendrocytes"))] <- "Glia"

# ==================================
# perform the actual deconvolution
# ==================================
bp_res_list <- list()

for (age in unique(bulk_meta$age)) {
  message( "processing: ", age)
  if (age == "P0") {
    ref_age <- "E18_S3"
  } else {
    ref_age <- age
  }
  
  i <- which( age_list == ref_age)

  fileName <- file_list[i]
  
  sc.data <- Read10X_h5(filename = paste( single_cell_base, "DiBella/", fileName, sep=""))
  
  sub_meta_data <- meta_data[ which( meta_data$biosample_id ==  age_list[i] & meta_data$New_cellType %in% cell_types_to_include ), ]
  
  # remove cell types with less than 100 cells
  cell_types_to_remove <- c( names ( which( table(sub_meta_data$New_cellType) < 50 ) ) )
  
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
  sub_meta_data$age <- ct_prefix[i]
  
  # ================================
  # deconvolution using BayesPrism
  # ================================
 
  sc.data <- as.matrix(sc.data)
  
  row_idx <- intersect( rownames(sc.data), gene_list )
  sc.data <- sc.data[row_idx, ]

  sc_dat_filtered <- cleanup.genes(
    input      = t(sc.data),
    input.type = "count.matrix",
    species    = "hs",            # "hs" for human, "mm" for mouse
    gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
    exp.cells  = 5                # exclude genes expressed in fewer than 5 cells
  )
  
  bulk_sample_list <- bulk_meta[ which( bulk_meta$age == age ), "sample_id"]
  
  # BayesPrism expects samples × genes
  bulk_matrix <- t(as.matrix( all[ , bulk_sample_list ] ))
  
  myPrism <- new.prism(
    reference         = sc_dat_filtered,
    mixture           = bulk_matrix,
    input.type        = "count.matrix",
    cell.type.labels  = sub_meta_data$broad_New_cellType,
    cell.state.labels = sub_meta_data$broad_New_cellType,
    key               = NULL,
    outlier.cut       = 0.01,
    outlier.fraction  = 0.1
  )
  
  bp_res_list[[ age ]] <- run.prism(
    prism   = myPrism,
    n.cores = 6        # parallelise — this step is slow
  )
  
}

# optional: save the object, since the processing takes a long time
saveRDS( object = bp_res_list, file = paste( main_base, "/RDS_files/decon_bp_res_list.rds", sep="" ) )
rownames( bulk_meta ) <- bulk_meta$sample_id
saveRDS( object = bulk_meta, file = paste( main_base, "/RDS_files/bulk_meta.rds", sep="" ) )

# ============================================
# start here to continue after deconvolution
# ============================================
bp_res_list <- readRDS( file = paste( main_base, "/RDS_files/decon_bp_res_list.rds", sep="" ) )
bulk_meta <- readRDS( file = paste( main_base, "/RDS_files/bulk_meta.rds", sep="" ) )

# =======================================================
# perform statistical analysis and prepare documentation
# =======================================================
list_of_plots <- list()

wb <- createWorkbook()

for (age in names( bp_res_list ) ) {
  
  props_bp <- get.fraction(
    bp        = bp_res_list[[ age ]],
    which.theta = "final",    # "final" = after convergence, "first" = initialisation
    state.or.type = "type"    # "type" for coarse, "state" for fine-grained
  )
  
  ct_diff_out <- run_coda_diffabundance (prop_mat = props_bp, group = bulk_meta[ rownames(props_bp), "genotype"], age = age)
  
  # write out raw proportions
  tmp_df <- as.data.frame( props_bp )
  tmp_df$genotype <- bulk_meta[ rownames(props_bp), "genotype"]
  
  addWorksheet(wb, paste( age, "BP_prop", sep="_" ))
  writeData( wb, paste( age, "BP_prop", sep="_" ), tmp_df )
  addFilter( wb, paste( age, "BP_prop", sep="_" ), row = 1, cols = 1:ncol( tmp_df ) )
  setColWidths( wb, paste( age, "BP_prop", sep="_" ), cols = 1:ncol( tmp_df ), widths="auto" )
  
  
  if ( length( ct_diff_out ) == 3 ) {
    
    tmp_df <- ct_diff_out[[3]]
    addWorksheet(wb, paste( age, "log_ratio", sep="_" ))
    writeData( wb, paste( age, "log_ratio", sep="_" ), tmp_df )
    addFilter( wb, paste( age, "log_ratio", sep="_" ), row = 1, cols = 1:ncol( tmp_df ) )
    setColWidths( wb, paste( age, "log_ratio", sep="_" ), cols = 1:ncol( tmp_df ), widths="auto" )
    
    list_of_plots[[ age ]] <- ct_diff_out[[ 2 ]]
    
    data.frame(
      statistic = ct_diff_out[[1]]$statistic,
      df        = ct_diff_out[[1]]$parameter,
      p.value   = ct_diff_out[[1]]$p.value,
      method    = ct_diff_out[[1]]$method
    ) -> tmp_df
    
    tmp_df$X1 <- rownames(tmp_df)
    
    addWorksheet(wb, age)
    writeData( wb, age, tmp_df )
    addFilter( wb, age, row = 1, cols = 1:ncol( tmp_df ) )
    setColWidths( wb, age, cols = 1:ncol( tmp_df ), widths="auto" )
    
  } else {
    
    list_of_plots[[ age ]] <- ct_diff_out[[ 5 ]]
    
    sheet_names <- paste(age, names( ct_diff_out[1:4] ), sep="_" )
  
    for (i in 1:4) {
      if (i == 1) {
        tmp_df <- as.data.frame( ct_diff_out[[ i ]]$stats )
      } else {
        tmp_df <- as.data.frame( ct_diff_out[[ i ]] )
      }
    
      tmp_df$X1 <- rownames(tmp_df)
      
      if (i==4) {
        tmp_df$genotype <- bulk_meta[ rownames(tmp_df), "genotype"]
      }
      addWorksheet(wb, sheet_names[i])
      writeData( wb, sheet_names[i], tmp_df )
      addFilter( wb, sheet_names[i], row = 1, cols = 1:ncol( tmp_df ) )
      setColWidths( wb, sheet_names[i], cols = 1:ncol( tmp_df ), widths="auto" )
      
    }
  }
}

saveWorkbook(wb, file = paste(main_base, "/Supplement/cellType_deconvolution_stats_1.xlsx", sep=""), overwrite = T) 

sessionInfo()
# R version 4.5.3 (2026-03-11)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS
# 
# Matrix products: default
# BLAS:   /opt/R/4.5.3/lib/R/lib/libRblas.so 
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.50.2               SummarizedExperiment_1.40.0 MatrixGenerics_1.22.0      
# [4] matrixStats_1.5.0           GenomicRanges_1.62.1        Seqinfo_1.0.0              
# [7] IRanges_2.44.0              S4Vectors_0.48.1            tibble_3.3.1               
# [10] dplyr_1.2.0                 tidyr_1.3.2                 limma_3.66.0               
# [13] compositions_2.0-9          zCompositions_1.6.2         truncnorm_1.0-9            
# [16] NADA_1.6-1.2                survival_3.8-6              MASS_7.3-65                
# [19] BayesPrism_2.2.3            NMF_0.28                    Biobase_2.70.0             
# [22] BiocGenerics_0.56.0         generics_0.1.4              cluster_2.1.8.2            
# [25] rngtools_1.5.2              registry_0.5-1              snowfall_1.84-6.3          
# [28] snow_0.4-4                  openxlsx_4.2.8.1            ggplot2_4.0.2              
# [31] data.table_1.18.2.1         Seurat_5.4.0                SeuratObject_5.3.0         
# [34] sp_2.2-1                   
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.23            splines_4.5.3               later_1.4.8                
# [4] bitops_1.0-9                polyclip_1.10-7             fastDummies_1.7.5          
# [7] lifecycle_1.0.5             edgeR_4.8.2                 doParallel_1.0.17          
# [10] globals_0.19.1              lattice_0.22-9              magrittr_2.0.4             
# [13] plotly_4.12.0               metapod_1.18.0              httpuv_1.6.17              
# [16] otel_0.2.0                  sctransform_0.4.3           spam_2.11-3                
# [19] zip_2.3.3                   spatstat.sparse_3.1-0       reticulate_1.45.0          
# [22] bayesm_3.1-7                cowplot_1.2.0               pbapply_1.7-4              
# [25] RColorBrewer_1.1-3          abind_1.4-8                 Rtsne_0.17                 
# [28] purrr_1.2.1                 tensorA_0.36.2.1            ggrepel_0.9.8              
# [31] irlba_2.3.7                 listenv_0.10.1              spatstat.utils_3.2-2       
# [34] goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.4-5      
# [37] dqrng_0.4.1                 fitdistrplus_1.2-6          parallelly_1.46.1          
# [40] codetools_0.2-20            DelayedArray_0.36.1         scuttle_1.20.0             
# [43] tidyselect_1.2.1            farver_2.1.2                ScaledMatrix_1.18.0        
# [46] spatstat.explore_3.8-0      jsonlite_2.0.0              BiocNeighbors_2.4.0        
# [49] progressr_0.18.0            ggridges_0.5.7              iterators_1.0.14           
# [52] foreach_1.5.2               tools_4.5.3                 ica_1.0-3                  
# [55] Rcpp_1.1.1                  glue_1.8.0                  gridExtra_2.3              
# [58] SparseArray_1.10.10         withr_3.0.2                 BiocManager_1.30.27        
# [61] fastmap_1.2.0               bluster_1.20.0              caTools_1.18.3             
# [64] digest_0.6.39               rsvd_1.0.5                  R6_2.6.1                   
# [67] mime_0.13                   colorspace_2.1-2            scattermore_1.2            
# [70] gtools_3.9.5                tensor_1.5.1                spatstat.data_3.1-9        
# [73] robustbase_0.99-7           httr_1.4.8                  htmlwidgets_1.6.4          
# [76] S4Arrays_1.10.1             uwot_0.2.4                  pkgconfig_2.0.3            
# [79] gtable_0.3.6                lmtest_0.9-40               S7_0.2.1                   
# [82] SingleCellExperiment_1.32.0 XVector_0.50.0              htmltools_0.5.9            
# [85] dotCall64_1.2               scales_1.4.0                png_0.1-9                  
# [88] spatstat.univar_3.1-7       scran_1.38.1                rstudioapi_0.18.0          
# [91] reshape2_1.4.5              nlme_3.1-168                zoo_1.8-15                 
# [94] stringr_1.6.0               KernSmooth_2.23-26          parallel_4.5.3             
# [97] miniUI_0.1.2                pillar_1.11.1               grid_4.5.3                 
# [100] vctrs_0.7.2                 RANN_2.6.2                  gplots_3.3.0               
# [103] promises_1.5.0              BiocSingular_1.26.1         beachmat_2.26.0            
# [106] xtable_1.8-8                packrat_0.9.3               cli_3.6.5                  
# [109] locfit_1.5-9.12             compiler_4.5.3              rlang_1.1.7                
# [112] future.apply_1.20.2         plyr_1.8.9                  stringi_1.8.7              
# [115] viridisLite_0.4.3           deldir_2.0-4                gridBase_0.4-7             
# [118] BiocParallel_1.44.0         lazyeval_0.2.2              spatstat.geom_3.7-3        
# [121] Matrix_1.7-5                RcppHNSW_0.6.0              patchwork_1.3.2            
# [124] future_1.70.0               statmod_1.5.2               shiny_1.13.0               
# [127] ROCR_1.0-12                 igraph_2.2.2                DEoptimR_1.2-0 


