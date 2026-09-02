# ==========================
# seventh step of analysis
# analysis of Ros Seq data
# =========================
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(plyr)

#for parallelisation of DESeq - use 7 cores
library(BiocParallel)
register(MulticoreParam(7))

set.seed(2401)

# =======================
# set up I/O folders
# =======================
base <- "./bulk_RNASeq/biorxiv/"

# ==============================================================
# prepare environment hash table as lookup for ENSMUSG - Symbol
# ==============================================================

#this comes from STAR index folder and was prepared solely by using info from Gencode
#initialise variables
ensmusg_symbol <- new.env()

#read conversion table
ensmusg_symbol_chr <- read.table(paste(base, "annotation/exchange_table_M27.tsv", sep=""), fill = T, header = F, stringsAsFactors = F)
colnames(ensmusg_symbol_chr) <- c("ENSMUSG", "SYMBOL", "chr")

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

# ==========================
# read the gene count data
# ==========================

# Folder with gene count files - give full path
gene_count_base=paste(base, "/Rptor_Ros_alignment/", sep="")

# EXEL file with meta data - give full path
fn_sample_list=paste(base, "/other_data/Ros_seq_sample_list.xlsx", sep="")

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
# add gene count file information to the meta data table
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
# Logs are availabel as Log.final.out.tar.gz on GitHub

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
ggplot( STAR_df2plot, aes(x=Condition, y=perc)) + geom_boxplot() + facet_grid(~ROX) + ylim(0,1)

ggsave( filename = paste(base, "/plots/QC/Ros_analysis.alignment_rate.pdf", sep="") )

# extract samples with unusually low alignment rate
idx2rem <- STAR_df2plot[which(as.numeric(sapply(STAR_df2plot$STAR.perc, function(x) gsub(pattern = "%", replacement = "", x=x))) < 50), "sample_id"]
# [1] "194050" "194052" "194054" "194058"
# prepare a group variable
sample_list_xls$group <- paste(sample_list_xls$Condition, sample_list_xls$ROX, sep="_")

# record in which step samples were removed:
round_rem <- list()
round_rem[[1]] <- idx2rem

rem <- c()
plot <- list()

# go through 2 rounds of QC and sample removal
for (round in 1:2) {
  
  sample_list_xls <- sample_list_xls[which(!sample_list_xls$sample_id %in% idx2rem),]
  
  # process all samples in comparisons - reduce noise by using samples with a mean read number < 10
  all_filter <- all[which(rowMeans(all[,2:ncol(all)]) > 10), ]
  all_filter <- all_filter[, c("gene_id", sapply(sample_list_xls$sample_id, as.character))]
  
  # prepare countdata - change class to integer (read counts have to be integer)
  countData <- sapply(all_filter[,c(2:(ncol(all_filter)))], as.integer)
  # add rownames - name of genes
  rownames(countData) <- all_filter[,c(1)]
  
  # prepare DESeq2 comparison - group by "sample" date in the design formula
  vsd <- varianceStabilizingTransformation(countData, blind=T)
  
  rv <- rowVars(vsd)
  ntop <- 500
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(vsd[select, ]))
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  d$sample <- rownames(d)
  
  d$Condition <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- paste(sample_list_xls$Condition[idx], sep="_")
    return(out)
  })
  
  d$ROX <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- sample_list_xls$ROX[idx]
    return(out)
  })
  
  plotParts <- ggplot(d, aes(x=PC1, y=PC2, color=Condition, shape=ROX)) + geom_point(size=3) + theme_classic()
  plotParts <- plotParts + ggtitle("PCA top 500 variable genes")
  plot[[length(plot)+1]] <- plotParts
  
  ggsave(filename = paste(base, "/plots/QC/Rox_analysis.", round, ".PCA_noCorrection_group.pdf", sep=""), plot = plot[[length(plot)]])
  
  # test for correlations between samples
  df2plot <- data.frame()
  for (group in unique(sample_list_xls$group)) {
    
    ids <- sample_list_xls[which(sample_list_xls$group == group), "sample_id"]
    pw_comp <- combn(ids, 2)
    
    for (i in 1: ncol(pw_comp)) {
      comp1 <- which(colnames(vsd) == pw_comp[1,i])
      comp2 <- which(colnames(vsd) == pw_comp[2,i])
      
      sample1 <- sample_list_xls[which(sample_list_xls$sample_id == pw_comp[1,i]), "sample_id"]
      sample2 <- sample_list_xls[which(sample_list_xls$sample_id == pw_comp[2,i]), "sample_id"]
      
      this_cor <- cor(vsd[,comp1], vsd[,comp2])
      
      df2plot <- rbind(df2plot, data.frame(group = group, cor=this_cor, sample1=sample1, sample2=sample2, stringsAsFactors = F))
    }
  }
  
  if (round == 1) {
    
    #PCA identifies one sample as an outlier
    idx2rem <- c(idx2rem, d[which(d$PC2 < -10), "sample"])
    
    #sample "194059" is an outlier
    round_rem[[2]] <- d[which(d$PC2 < -10), "sample"]
    
 }
  
  plot[[length(plot)+1]] <- ggplot(df2plot, aes(x=group, y=cor, group=group)) + geom_boxplot() + geom_point() + theme_classic() + ylim(0,1) + geom_hline(yintercept = 0.9)
  ggsave(filename = paste(base, "/plots/QC/Rox_correlation." , round, ".corr_pw_comp_samples.pdf", sep=""), plot = plot[[length(plot)]])
  
}

sample_list_xls <- sample_list_xls[which(!sample_list_xls$sample_id %in% idx2rem),]
sample_list_xls$Condition <- gsub(pattern = " ", replacement = "_", x = sample_list_xls$Condition)

# sample numbers:
ddply(sample_list_xls, c("ROX", "Condition"), summarise,
               N    = length(Condition)
)

#   ROX    Condition N
# 1 neg Control-MADM 5
# 2 neg   Rptor-MADM 4
# 3 pos Control-MADM 3
# 4 pos   Rptor-MADM 3

# ========================================================================
# Note:
# Rptor is too lowly expressed to be detected by this shallow sequencing
# no check for deletion allele possible
# ========================================================================

# ========================================
# DEG analysis
# ========================================

res <- list()
norm_counts <- list()
vsd_list <- list()

for( rox in unique(sample_list_xls$ROX)) {

    sl_txt_sub <- sample_list_xls[which(sample_list_xls$ROX == rox),]
                                        
    # extract the columns from expression data matching the list of samples extracted above
    # also ectract gene_id!
    colnames_to_extract <- c("gene_id", sl_txt_sub$sample_id)
    
    # NOTE: expression cutoff used here is 10
    all_filter <- all[, which(colnames(all) %in% colnames_to_extract)]
    all_filter <- all_filter[which(rowMeans(all_filter[,2:ncol(all_filter)]) > 10),]
    
    # prepare metadata for later DESeq2 analysis    
    meta<-data.frame(  cond=sl_txt_sub$Condition, rox=sl_txt_sub$ROX )
    rownames(meta) <- colnames(all_filter[,2:(ncol(all_filter))])
    
    # prepare countdata - change class to integer (read counts have to be integer)
    countData <- sapply(all_filter[,c(2:(ncol(all_filter)))], as.integer)
    # add rownames - name of genes
    rownames(countData) <- all_filter[,c(1)]
    
    # prepare DESeq2 comparison
    ddsMat<-DESeqDataSetFromMatrix(countData = countData,
                                   colData=meta,
                                   design= ~ cond)
    
    # do DESeq2 analysis
    dds <- DESeq(ddsMat, fitType="local", quiet=T, betaPrior = F, parallel = T)

    # save the normalised counts
    norm_counts[[rox]] <- counts(dds, normalize=T)
    
    # prepare variance stabilised counts
    vsd <- vst(dds, blind=FALSE)
    vsd_list[[rox]] <- vsd
    
    # prepare a list of contrasts to be performed
    groups_contrasts <- list()
    groups_contrasts[[1]] <- c("cond", "Rptor-MADM", "Control-MADM")

    for (x in 1:length(groups_contrasts)) {
      cat(x, "\n")
      descr <- paste( groups_contrasts[[x]][2], groups_contrasts[[x]][3], sep="_")
      this_res <- as.data.frame(results(dds, contrast = groups_contrasts[[x]]))
      this_res <- this_res[order(this_res$padj, decreasing = F),]
      this_res <- this_res[order(this_res$padj),]
      this_res$symbol <- sapply(rownames(this_res), function (x) ensmusg_symbol[[x]])
      res[[paste(rox, descr, sep="_")]] <- this_res
    }
}

# ================================
# cell type enrichment analysis
# ================================

# read cell type marker genes
CT_GENE_DF <- readRDS( file = paste(base, "/other_data/DiBella.P1.marker_df.RDS", sep="") )

# combine gsea analysis for plotting
gsea_res_df <- data.frame()
# save gsea analysis for reporting
gsea_res_list <- list() 
for ( comp in c("pos_Rptor-MADM_Control-MADM", "neg_Rptor-MADM_Control-MADM")) {
  
  # obtain DEG analysis and remove duplicated symbols
  tmp_de <- res[[ comp ]]
  tmp_de <- tmp_de[!(duplicated(tmp_de$symbol) | duplicated(tmp_de$symbol, fromLast = TRUE)), ]
  
  # create a gene list ordered by a score that is defined by log10( pvalue )
  # prefix corrected for fold change
  tmp_de <- tmp_de[ which(!is.na(tmp_de$pvalue)), ]
  tmp_de$score <- log10(tmp_de$pvalue) * -1
  tmp_de$score[ which(tmp_de$log2FoldChange > 0) ] <- tmp_de$score[ which(tmp_de$log2FoldChange > 0) ] * -1 
  tmp_de <- tmp_de[ order(tmp_de$score, decreasing = T), ]
  
  de_vec <- tmp_de$score
  names(de_vec) <- tmp_de$symbol
  
  gsea_result <- GSEA(de_vec, TERM2GENE = CT_GENE_DF, seed = T, pvalueCutoff = 0.9)
  gsea_res_list[[ comp ]] <- gsea_result@result
  
  # extract terms and combine data in a way I can plot a heatmap later
  gsea_result <- gsea_result@result[, c(1,5,7)]
  gsea_result$score <- log10(gsea_result$p.adjust)
  gsea_result$score[which( gsea_result$NES > 0)] <- gsea_result$score[which( gsea_result$NES > 0)] * -1
  
  colnames(gsea_result) <- c("ID", "NES", "padj", comp)
  
  if (ncol(gsea_res_df) == 0){
    gsea_res_df <- gsea_result[,c(1,4)]
  } else {
    gsea_res_df <- merge(gsea_res_df, gsea_result[,c(1,4)], by="ID", all=T)
  }
}

# very low number of ties - accept for this analysis

# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# preparing geneSet collections...
# GSEA analysis...
# leading edge analysis...
# done...
# Warning message:
#   In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                There are ties in the preranked stats (0.01% of the list).
#                              The order of those tied genes will be arbitrary, which may produce unexpected results.
                             
# convert df to a matrix, cut values for easier plotting
rownames( gsea_res_df ) <- gsea_res_df$ID
mat_score <- gsea_res_df[ ,c(2:ncol(gsea_res_df)) ]
mat_score[ is.na(mat_score) ] <- 0

# create a plot only for UL and DL neurons
df2plot <- reshape2::melt( mat_score[c("UL_CPN_P1", "DL_CPN_P1"),] )
df2plot$cell_type <- rep( rownames(mat_score[c("UL_CPN_P1", "DL_CPN_P1"),]), ncol(mat_score[c("UL_CPN_P1", "DL_CPN_P1"),]))
df2plot$cell_type <- factor( df2plot$cell_type, levels = c("UL_CPN_P1", "DL_CPN_P1"))

ggplot( df2plot, aes(x=variable, y=abs(value), fill=cell_type)) + geom_bar(stat="identity", position="dodge", color="black") + 
  theme_classic() + ylab("significance score") + xlab("") + 
  scale_fill_manual( values=c("DL_CPN_P1" = "grey20", "UL_CPN_P1" = "white")) +
  ggtitle("cell type marker enrichment\nRptor-MADM/Rptor-Control") +
  scale_x_discrete( labels = c("Ros positive", "Ros negative")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,3) +
  geom_hline(yintercept = 1.3) # significance cutoff - equals to adjusted  p-value of 0.05

ggsave(filename = paste(base, "/plots/Supplement/FigS9D_Rox_cellType_enrichment.pdf", sep=""), width=3)

# =================
# final reporting
# =================

# sample list
sample_list_xls$used_for_analysis <- "Y"
sample_list_reporting <- merge( STAR_df2plot, sample_list_xls[,c("sample_id", "used_for_analysis")], by="sample_id", all=T)
sample_list_reporting$used_for_analysis[ which( is.na(sample_list_reporting$used_for_analysis) ) ] <- "N"
sample_list_reporting <- sample_list_reporting[,c(1:8,11)]

wb <- createWorkbook()

addWorksheet(wb, "Rox sample list")
writeData(wb, "Rox sample list", sample_list_reporting)
addFilter(wb, "Rox sample list", row = 1, cols = 1:ncol( sample_list_reporting ))
setColWidths(wb, "Rox sample list", cols = 1:ncol(sample_list_reporting), widths="auto")

saveWorkbook(wb, file = paste(base, "Supplement/Rox_sample_list.xlsx", sep=""), overwrite = T) 

# DEG analysis
wb <- createWorkbook()

for (name in names(res)) {

  addWorksheet(wb, name)
  
  writeData(wb, name, res[[name]])
  addFilter(wb, name, row = 1, cols = 1:ncol(res[[name]]))
  setColWidths(wb, name, cols = 1:ncol(res[[name]]), widths="auto")
}

saveWorkbook(wb, file = paste(base, "Supplement/Rox_DEG.xlsx", sep=""), overwrite = T) 

# GSEA analysis

wb <- createWorkbook()

for (name in names(gsea_res_list)) {
  
  addWorksheet(wb, name)
  writeData(wb, name, gsea_res_list[[name]])
  addFilter(wb, name, row = 1, cols = 1:ncol(gsea_res_list[[name]]))
  setColWidths(wb, name, cols = 1:ncol(gsea_res_list[[name]]), widths="auto")
}

saveWorkbook(wb, file = paste(base, "Supplement/Rox_GSEA.xlsx", sep=""), overwrite = T) 

# data for the actual plot in the figure
wb <- createWorkbook()

addWorksheet(wb, "Rox GSEA plot")
writeData(wb, "Rox GSEA plot", df2plot)
addFilter(wb, "Rox GSEA plot", row = 1, cols = 1:ncol( df2plot ))
setColWidths(wb, "Rox GSEA plot", cols = 1:ncol(df2plot), widths="auto")

saveWorkbook(wb, file = paste(base, "Supplement/Rox_gsea_bar_raw.xlsx", sep=""), overwrite = T) 

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
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] BiocParallel_1.36.0         plyr_1.8.9                  clusterProfiler_4.10.1      ggplot2_3.5.0              
# [5] DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0              MatrixGenerics_1.14.0      
# [9] matrixStats_1.2.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.8         IRanges_2.36.0             
# [13] S4Vectors_0.40.2            BiocGenerics_0.48.1         openxlsx_4.2.5.2           
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.2               bitops_1.0-7            gson_0.1.0              shadowtext_0.1.3        gridExtra_2.3          
# [6] rlang_1.1.3             magrittr_2.0.3          DOSE_3.28.2             compiler_4.3.2          RSQLite_2.3.6          
# [11] systemfonts_1.0.6       png_0.1-8               vctrs_0.6.5             reshape2_1.4.4          stringr_1.5.1          
# [16] pkgconfig_2.0.3         crayon_1.5.2            fastmap_1.1.1           XVector_0.42.0          labeling_0.4.3         
# [21] ggraph_2.2.1            utf8_1.2.4              HDO.db_0.99.1           enrichplot_1.22.0       ragg_1.3.0             
# [26] purrr_1.0.2             bit_4.0.5               zlibbioc_1.48.2         cachem_1.0.8            aplot_0.2.2            
# [31] jsonlite_1.8.8          blob_1.2.4              DelayedArray_0.28.0     tweenr_2.0.3            parallel_4.3.2         
# [36] R6_2.5.1                RColorBrewer_1.1-3      stringi_1.8.3           GOSemSim_2.28.1         Rcpp_1.0.12            
# [41] Matrix_1.6-5            splines_4.3.2           igraph_2.1.4            tidyselect_1.2.1        viridis_0.6.5          
# [46] qvalue_2.34.0           rstudioapi_0.16.0       abind_1.4-5             codetools_0.2-19        lattice_0.21-9         
# [51] tibble_3.2.1            treeio_1.26.0           withr_3.0.0             KEGGREST_1.42.0         gridGraphics_0.5-1     
# [56] scatterpie_0.2.2        polyclip_1.10-6         zip_2.3.1               Biostrings_2.70.3       ggtree_3.10.1          
# [61] pillar_1.9.0            packrat_0.9.2           ggfun_0.1.4             generics_0.1.3          RCurl_1.98-1.14        
# [66] tidytree_0.4.6          munsell_0.5.0           scales_1.3.0            glue_1.7.0              lazyeval_0.2.2         
# [71] tools_4.3.2             data.table_1.15.0       fgsea_1.28.0            locfit_1.5-9.12         fs_1.6.3               
# [76] graphlayouts_1.1.1      fastmatch_1.1-4         tidygraph_1.3.1         cowplot_1.1.3           grid_4.3.2             
# [81] ape_5.8                 tidyr_1.3.1             AnnotationDbi_1.64.1    colorspace_2.1-0        nlme_3.1-163           
# [86] patchwork_1.2.0         GenomeInfoDbData_1.2.11 ggforce_0.4.2           cli_3.6.2               textshaping_0.3.7      
# [91] fansi_1.0.6             viridisLite_0.4.2       S4Arrays_1.2.1          dplyr_1.1.4             gtable_0.3.4           
# [96] yulab.utils_0.1.4       digest_0.6.34           ggplotify_0.1.2         SparseArray_1.2.4       ggrepel_0.9.5          
# [101] farver_2.1.1            memoise_2.0.1           lifecycle_1.0.4         httr_1.4.7              GO.db_3.18.0           
# [106] bit64_4.0.5             MASS_7.3-60 
