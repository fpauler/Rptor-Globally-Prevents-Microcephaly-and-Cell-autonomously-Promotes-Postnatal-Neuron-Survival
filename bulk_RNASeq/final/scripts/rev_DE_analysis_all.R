# ===================================================
# analyse bulk RNA-Seq without sample removal
# this follows a reviewer's question:
# originally one P0 WT sample was removed from analysis
# test here what impact that removal had on the analysis
# ===================================================

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (clusterProfiler)
library (org.Mm.eg.db)
library (pheatmap)
library (dplyr)
library (ggplot2)
library (tibble)
library (tidyr)


# for parallelisation of DESeq - use 7 cores
library(BiocParallel)
register(MulticoreParam(7))

set.seed(2401)

bulk_base <- "./bulk_RNASeq/biorxiv/" # initial analysis
main_base <- "./bulk_RNASeq/final/" # main working folder

# ==============================================================
# prepare environment hash table as lookup for ENSMUSG - Symbol
# ==============================================================

# this comes from STAR index folder and was prepared solely by using info from Gencode
# initialise variables
ensmusg_symbol <- new.env()
ensmusg_chr <- new.env()
ensmusg_entrezid <- new.env()

# read conversion table - this part of the analysis was performed with:
# genome sequence (GRCm38.p5), gene annotation conversion tables and Gencode M16

ensmusg_symbol_chr <- read.table(paste(bulk_base, "annotation/exchange.table_M16.tsv", sep=""), fill = T, header = F, stringsAsFactors = F)
colnames(ensmusg_symbol_chr) <- c("ENSMUSG", "ENSMUST", "SYMBOL", "chr", "ENTREZID")

# fill conversion table - handled just like a list but apparently much faster
for (x in 1:nrow(ensmusg_symbol_chr)) {
  # deal here with multi mappers - one SYMBOL mapped to many ENSMUSG
  # a comma in the SYMBOL name indicates these multi-mappers
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
fn_sample_list=paste(bulk_base, "/other_data/Rptor_seq_sample_list.xlsx", sep="")

# read EXEL file using openxlsx
sample_list_xls <- read.xlsx(xlsxFile = fn_sample_list, sheet = 1, skipEmptyRows = FALSE)

# match samples with read count files - done via sampleID
# recursive search for all files in the count folder - returns relative path
files <- list.files(path = gene_count_base, pattern="*ReadsPerGene.out.tab", recursive = T)
# make path absolute
full_fn_path_list <- lapply(files, function (x) {paste(gene_count_base, x, sep="")})

# prepare dataframe from file list - matching ID with filename
files.DF <- data.frame(fn = unlist(full_fn_path_list), stringsAsFactors = F)
files.DF$sample_id <- apply(files.DF, 1, function (x) {
  # extract samle ID from filename
  parts <- strsplit(x[1], "\\.")
  return(parts[[1]][length(parts[[1]])-3])
})

sample_list_xls <- merge(sample_list_xls, files.DF, by="sample_id")

# initialise dataframe holding all read counts df
all <- data.frame(stringsAsFactors = F)

#read in all counts and merge in one dataframe - by gene_id
for (line in 1:nrow(sample_list_xls)) {
  this_count <- read.table(file = sample_list_xls$fn[line], header = F, sep="\t", stringsAsFactors = F)
  # remove first 4 lines - only basic stats of features
  this_count<-this_count[5:nrow(this_count),c(1,2)]
  
  colnames(this_count) <- c("gene_id", line)
  
  if (nrow(all) == 0) {
    all <- this_count
  } else {
    all <- merge(all, this_count, by="gene_id")
  }
}

# give meaningful column names - sample_id in this case
colnames(all)<-c("gene_id", sample_list_xls$sample_id)

# ========
# QC
# ========

# analyse STAR Logs
files.STAR <- c(list.files(path=gene_count_base, pattern=glob2rx("*.Log.final.out"), recursive = T))

STAR_df2plot <- data.frame()
for (file_name in files.STAR) {
  STAR.Log <- read.table(file = paste(gene_count_base, file_name, sep=""), header = F, fill = T, sep="\t", stringsAsFactors = F)
  
  # isolate alignment folder name
  folders <- strsplit(x = file_name, split = "\\/")[[1]] 
  
  # if the folder name contains a dot this indicates that the sample ID needs to be extracted
  # otherwise the sample_id is in the filename
  sample_name <- strsplit(x = folders[[length(folders)-1]], split = "\\.")[[1]][2]

  STAR_df2plot <- rbind(STAR_df2plot, data.frame(sample_id=sample_name, STAR.unique.aligned=STAR.Log$V2[8], 
                                                 STAR.tot=STAR.Log$V2[5], STAR.perc = STAR.Log$V2[9], stringsAsFactors = F))
}

STAR_df2plot <- merge(STAR_df2plot, sample_list_xls, by="sample_id")
STAR_df2plot$perc <- as.numeric( STAR_df2plot$STAR.unique.aligned ) / as.numeric( STAR_df2plot$STAR.tot )

# identical to the first run - just for double check
ggplot( STAR_df2plot, aes(x=genotype, y=perc)) + geom_boxplot() + facet_grid(~age) + ylim(0,1)

# no samples with unusually low alignment rate, as before
STAR_df2plot[which(as.numeric(sapply(STAR_df2plot$STAR.perc, function(x) gsub(pattern = "%", replacement = "", x=x))) < 50), "sample_name"]
# NULL

# prepare a group variable
sample_list_xls$group <- paste(sample_list_xls$age, sample_list_xls$genotype, sep="_")

# sample numbers, only difference is the P0 WT sample number
sample_list_xls %>% group_by( age, genotype ) %>% summarise ( n=n() )

# `summarise()` has grouped output by 'age'. You can override using the `.groups` argument.
# # A tibble: 9 × 3
# # Groups:   age [3]
# age   genotype     n
# <chr> <chr>    <int>
# 1 E13   HET          4
# 2 E13   WT           7
# 3 E13   cKO          4
# 4 E16   HET          4
# 5 E16   WT           5
# 6 E16   cKO          4
# 7 P0    HET          4
# 8 P0    WT           3 <- was 2 before!
# 9 P0    cKO          5

files <- list.files(path = gene_count_base, pattern="Rptor_exon6.counts.out", recursive = T)
# make path absolute
full_fn_path_list <- sapply(files, function (x) {paste(gene_count_base, x, sep="")})

delcounts <- data.frame()
for (file in files) {
  seq_id <- strsplit(x = file, split = "/|\\.")[[1]][2]
  if (seq_id %in% sample_list_xls$sample_id) {
    tmp <- read.table(file = paste(gene_count_base, file, sep=""), header = F, stringsAsFactors = F)
    tmp <- tmp[,c("V6", "V7")]
    colnames(tmp) <- c("gene_id", seq_id)
    if (nrow(delcounts) == 0){
      delcounts <- tmp
    } else {
      delcounts <- merge(delcounts, tmp, by="gene_id")
    }
  }
}

delcounts$gene_id[which(delcounts$gene_id == "Rptor_ex6")] <- "Rptor_del"

res <- list()
norm_counts <- list()
vsd_list <- list()

# ========================================
# DEG analysis - separately for each age
# ========================================

# DEG analysis - separately for each age
for(age in c("E13", "E16", "P0")) {

    # extracts samples to analyse - based an age
    sl_txt_sub <- sample_list_xls[which(sample_list_xls$age == age),]
                                        
    # extract the columns from expression data fitting to the list of samples extracted above
    # also extract gene_id!
    colnames_to_extract <- c("gene_id", sl_txt_sub$sample_id)
    
    # NOTE: expression cutoff used here is 10
    all_filter <- all[, which(colnames(all) %in% colnames_to_extract)]
    all_filter <- all_filter[which(rowMeans(all_filter[,2:ncol(all_filter)]) > 10),]
    
    # add delcounts data
    all_filter <- rbind(all_filter, delcounts[,colnames(all_filter)])
    
    # prepare metadata for later DESeq2 analysis    
    meta<-data.frame(  gene=sl_txt_sub$gene, genotype=sl_txt_sub$genotype, sample = sl_txt_sub$sample, age = sl_txt_sub$age, color=sl_txt_sub$color)
    rownames(meta) <- colnames(all_filter[,2:(ncol(all_filter))])
    
    # prepare countdata - change class to integer (read counts have to be integer)
    countData <- sapply(all_filter[,c(2:(ncol(all_filter)))], as.integer)
    # add rownames - name of genes
    rownames(countData) <- all_filter[,c(1)]
    
    # prepare DESeq2 comparison
    ddsMat<-DESeqDataSetFromMatrix(countData = countData,
                                   colData=meta,
                                   design= ~ genotype)
    
    # do DESeq2 analysis
    dds <- DESeq(ddsMat, fitType="local", quiet=T, betaPrior = F, parallel = T)

    # save the normalised counts
    norm_counts[[age]] <- counts(dds, normalize=T)
    
    # prepare variance stabilised counts
    vsd <- vst(dds, blind=FALSE)
    vsd_list[[age]] <- vsd
    
    # prepare a list of contrasts to be performed
    groups_contrasts <- list()
    groups_contrasts[[1]] <- c("genotype", "HET", "WT")
    groups_contrasts[[2]] <- c("genotype", "cKO", "WT")
    groups_contrasts[[3]] <- c("genotype", "cKO", "HET")
    
    # comparisons of genotypes using contrasts
    for (x in 1:length(groups_contrasts)) {
      cat(x, "\n")
      descr <- paste( groups_contrasts[[x]][2], groups_contrasts[[x]][3], sep="_")
      this_res <- as.data.frame(results(dds, contrast = groups_contrasts[[x]]))
      this_res <- this_res[order(this_res$padj, decreasing = F),]
      this_res <- this_res[order(this_res$padj),]
      this_res$symbol <- sapply(rownames(this_res), function (x) ensmusg_symbol[[x]])
      res[[paste(age, descr, sep="_")]] <- this_res
    }
}

# GO analysis
# load filtered GO term df
GO_GENE_DF <- readRDS( file = paste(bulk_base, "/other_data/GO_Term_df.RDS", sep="") )
GO_GENE_DF <- GO_GENE_DF[,c(2,1)]

# basic DEG stats
eGO <- list()
deg_list <- list()
deg_count <- list()

# ==================
# compare analyses
# this step prepares comparison of the main findings between:
# old: one WT sample removed - analysis as is in the manuscript
# new: no sample removal
# ==================

# read the DEG analysis, filter DEGs and do GO term enrichment analysis

age <- "P0"
for (analysis in c("old", "new")) {
  for ( group in c("cKO", "HET") ) {
    if (analysis == "old") {
      deg_df <- read.xlsx( paste(bulk_base, "Supplement/DEG_raw_DEG_stats.xlsx", sep=""), sheet = paste( age, group, "WT", sep="_" ))
    } else {
      deg_df <- res[[  paste( age, group, "WT", sep="_" ) ]]
    }
    
    deg_list[[ paste( analysis, age, group, "all", sep="_") ]] <- deg_df$symbol[ which(deg_df$padj < 0.05) ]
    
    eGO[[ paste( analysis, age, group, sep="_") ]] <- NA
    eGO[[ paste( analysis, age, group, sep="_") ]] <-  enricher ( gene = deg_list[[ paste( analysis, age, group, "all", sep="_") ]], 
                                                        TERM2GENE = GO_GENE_DF, pvalueCutoff = 0.9, qvalueCutoff = 0.9,
                                                        minGSSize = 10, 
                                                        maxGSSize = 2000
    )
  }
}

# ================================================================================
# determine overlap of Rptor-MADM/Rptor-Control and cKO-Rptor-MADM/Rptor-Control
# in both old and new analysis
# ================================================================================

library(ggVennDiagram)
library (cowplot)

plot_list <- list()

for (analysis in c("old", "new")) {
  
  tmp_deg_list <- list(c(deg_list[[ paste(analysis, age, "cKO_all", sep="_") ]]), c(deg_list[[ paste(analysis, age, "HET_all", sep="_") ]]) )
  names(tmp_deg_list) <- c(paste(age, "cKO", sep="_"), paste(age, "het", sep="_") )

  plot_list[[ analysis ]] <- ggVennDiagram( x = tmp_deg_list, label_geom = "label", label = "count", label_alpha = 0 ) + 
    scale_fill_gradient(low="white",high = "white") +
    theme(legend.position = "NONE") + coord_cartesian()
}

venn_plot <- plot_grid( plotlist = plot_list, ncol = 1, greedy = F )
ggsave  ( filename = paste(main_base, "/plots/reviewer/Rev_Fig_DEG_venn_plot_comp.pdf", sep=""), plot = venn_plot, width=3, height=5)

# ==========================
# comparison of GO analyses
# ==========================

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

GOgroup_df2plot <- data.frame()

# note that both analyses (old and new) are stored in eGO
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
    
    GOgroup_df2plot <- rbind( GOgroup_df2plot, data.frame( GO_group = group, value = max_score, analysis = id_parts[1], age = id_parts[2], genotype = id_parts[3]) )
  }
  
}

GOgroup_df2plot$value[which(GOgroup_df2plot$value > 2.5)] <- 2.5 
GOgroup_df2plot$GO_group <- factor( GOgroup_df2plot$GO_group, levels = rev( c("cell_cycle", "cytokine", "TOR", "reactive_oxygen", "metabol", "mitoch", "axon", "apop")))
GOgroup_df2plot$genotype <- factor( GOgroup_df2plot$genotype, levels = c("HET", "cKO"))

ggplot(GOgroup_df2plot, aes(x=genotype, y=GO_group, fill=value)) + geom_tile() + 
  scale_fill_gradient(low="white", high = "black", limits = c(0,2.5)) + 
  facet_wrap(~analysis) + theme_classic() + ggtitle("transcriptomic GO term analysis")
ggsave  ( filename = paste(main_base, "/plots/reviewer/Rev_Fig_DEG_GO_score_plot.pdf", sep=""), width=4, height=4)

# ============================================================
# plot the number of DEGs in the "old" and the "new" analysis
# ============================================================

df2plot <- as.data.frame( sapply(deg_list, length) )
colnames(df2plot) <- "nDEG"

df2 <- df2plot |>
  rownames_to_column("row_id") |>
  separate(
    row_id,
    into = c("analysis", "age", "genotype", "type"),  # change to your metadata names
    sep = "_",
    remove = FALSE
  )

ggplot( df2, aes(x=genotype, y=nDEG)) + geom_bar(stat="identity", position="dodge") + 
  facet_grid(~analysis) + theme_classic()
ggsave  ( filename = paste(main_base, "/plots/reviewer/Rev_Fig_number_DEG_plot.pdf", sep=""), width=4, height=4)

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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8        LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8   
# [6] LC_MESSAGES=en_GB.UTF-8    LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DESeq2_1.50.2               SummarizedExperiment_1.40.0 MatrixGenerics_1.22.0       matrixStats_1.5.0           GenomicRanges_1.62.1       
# [6] Seqinfo_1.0.0               ggplot2_4.0.2               BiocParallel_1.44.0         pheatmap_1.0.13             org.Mm.eg.db_3.22.0        
# [11] AnnotationDbi_1.72.0        IRanges_2.44.0              S4Vectors_0.48.1            Biobase_2.70.0              BiocGenerics_0.56.0        
# [16] generics_0.1.4              clusterProfiler_4.18.4      openxlsx_4.2.8.1           
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.3.0               gson_0.1.0              rlang_1.1.7             magrittr_2.0.4          DOSE_4.4.0              otel_0.2.0             
# [7] compiler_4.5.3          RSQLite_2.4.6           png_0.1-9               systemfonts_1.3.2       vctrs_0.7.2             reshape2_1.4.5         
# [13] stringr_1.6.0           pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0           XVector_0.50.0          labeling_0.4.3         
# [19] utf8_1.2.6              enrichplot_1.30.5       ragg_1.5.2              purrr_1.2.1             bit_4.6.0               cachem_1.1.0           
# [25] aplot_0.2.9             jsonlite_2.0.0          blob_1.3.0              tidydr_0.0.6            DelayedArray_0.36.1     tweenr_2.0.3           
# [31] cluster_2.1.8.2         parallel_4.5.3          R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3      GOSemSim_2.36.0        
# [37] Rcpp_1.1.1              ggtangle_0.1.2          R.utils_2.13.0          Matrix_1.7-5            splines_4.5.3           igraph_2.2.2           
# [43] tidyselect_1.2.1        abind_1.4-8             qvalue_2.42.0           rstudioapi_0.18.0       codetools_0.2-20        lattice_0.22-9         
# [49] tibble_3.3.1            plyr_1.8.9              withr_3.0.2             treeio_1.34.0           KEGGREST_1.50.0         S7_0.2.1               
# [55] gridGraphics_0.5-1      scatterpie_0.2.6        polyclip_1.10-7         zip_2.3.3               Biostrings_2.78.0       pillar_1.11.1          
# [61] BiocManager_1.30.27     ggtree_4.0.5            packrat_0.9.3           ggfun_0.2.0             scales_1.4.0            tidytree_0.4.7         
# [67] glue_1.8.0              gdtools_0.5.0           lazyeval_0.2.2          tools_4.5.3             ggnewscale_0.5.2        data.table_1.18.2.1    
# [73] locfit_1.5-9.12         fgsea_1.36.2            ggiraph_0.9.6           fs_2.0.1                fastmatch_1.1-8         cowplot_1.2.0          
# [79] grid_4.5.3              tidyr_1.3.2             ape_5.8-1               nlme_3.1-168            patchwork_1.3.2         ggforce_0.5.0          
# [85] cli_3.6.5               rappdirs_0.3.4          textshaping_1.0.5       fontBitstreamVera_0.1.1 S4Arrays_1.10.1         dplyr_1.2.0            
# [91] gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.2.4       digest_0.6.39           fontquiver_0.2.1        SparseArray_1.10.10    
# [97] ggrepel_0.9.8           ggplotify_0.1.3         htmlwidgets_1.6.4       farver_2.1.2            memoise_2.0.1           htmltools_0.5.9        
# [103] R.oo_1.27.1             lifecycle_1.0.5         httr_1.4.8              GO.db_3.22.0            fontLiberation_0.1.0    bit64_4.6.0-1          
# [109] MASS_7.3-65  
