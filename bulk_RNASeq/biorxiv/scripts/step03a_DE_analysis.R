# ===================================================
# third step of the analysis
# analyse bulk RNA-Seq up to the stage of DEG stats
# ===================================================

library(openxlsx)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(pheatmap)
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
ensmusg_chr <- new.env()
ensmusg_entrezid <- new.env()

# read conversion table - this part of the analysis was performed with:
# genome sequence (GRCm38.p5), gene annotation conversion tables and Gencode M16

ensmusg_symbol_chr <- read.table(paste(base, "annotation/exchange.table_M16.tsv", sep=""), fill = T, header = F, stringsAsFactors = F)
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
gene_count_base=paste(base, "/Rptor_SMARTer.files/", sep="")

# EXEL file with meta data - give full path
fn_sample_list=paste(base, "/other_data/Rptor_seq_sample_list.xlsx", sep="")

# read EXEL file using openxlsx
sample_list_xls <- read.xlsx(xlsxFile = fn_sample_list, sheet = 1, skipEmptyRows = FALSE)

# match meta data with read count files - done via sampleID
# recursive search for all files in the count folder - returns relative path
files <- list.files(path = gene_count_base, pattern="*ReadsPerGene.out.tab", recursive = T)
# make path absolute
full_fn_path_list <- lapply(files, function (x) {paste(gene_count_base, x, sep="")})
# prepare dataframe from file list - matching ID with filename
files.DF <- data.frame(fn = unlist(full_fn_path_list), stringsAsFactors = F)
files.DF$sample_id <- apply(files.DF, 1, function (x) {
  # extract sample ID from filename
  parts <- strsplit(x[1], "\\.")
  return(parts[[1]][length(parts[[1]])-3])
})

sample_list_xls <- merge(sample_list_xls, files.DF, by="sample_id")

# initialise dataframe
all <- data.frame(stringsAsFactors = F)

#read in all counts and merge in one dataframe - by gene_id
for (line in 1:nrow(sample_list_xls)) {
  this_count <- read.table(file = sample_list_xls$fn[line], header = F, sep="\t", stringsAsFactors = F)
  
  #remove first 4 lines - only basic stats of features
  this_count<-this_count[5:nrow(this_count),c(1,2)]
  
  colnames(this_count) <- c("gene_id", line)
  
  if (nrow(all) == 0) {
    all <- this_count
  } else {
    all <- merge(all, this_count, by="gene_id")
  }
}

# give meaningful column names - sample_id in this case
colnames(all) <- c("gene_id", sample_list_xls$sample_id)

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

ggplot( STAR_df2plot, aes(x=genotype, y=perc)) + geom_boxplot() + facet_grid(~age) + ylim(0,1)

# no apparent difference between the genotypes
ggsave( filename = paste(base, "/plots/QC/DEG_analysis.alignment_rate.pdf", sep="") )

# extract samples with unusually low alignment rate
idx2rem <- STAR_df2plot[which(as.numeric(sapply(STAR_df2plot$STAR.perc, function(x) gsub(pattern = "%", replacement = "", x=x))) < 50), "sample_name"]

# prepare a group variable
sample_list_xls$group <- paste(sample_list_xls$age, sample_list_xls$genotype, sep="_")

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
  
  vsd <- varianceStabilizingTransformation(countData, blind=T)
  
  rv <- rowVars(vsd)
  ntop <- 500
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(vsd[select, ]))
  
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  d$sample <- rownames(d)
  
  d$age <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- paste(sample_list_xls$age[idx], sep="_")
    return(out)
  })
  
  d$genotype <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- sample_list_xls$genotype[idx]
    return(out)
  })
  
  d$date <- sapply(d$sample, function (x) {
    idx <- which(sample_list_xls$sample_id == x)
    out <- sample_list_xls$sample[idx]
    return(out)
  })
  
  plotParts <- ggplot(d, aes(x=PC1, y=PC2, color=age, shape=genotype)) + geom_point(size=3) + theme_classic()
  plotParts <- plotParts + ggtitle("PCA top 500 variable genes")
  plot[[length(plot)+1]] <- plotParts
  
  ggsave(filename = paste(base, "/plots/QC/analysis.", round, ".PCA_noCorrection_group.pdf", sep=""), plot = plot[[length(plot)]])
  
  plotParts <- ggplot(d, aes(x=PC1, y=PC2, color=as.factor(date), shape=genotype)) + geom_point(size=3) + theme_classic()
  plotParts <- plotParts + ggtitle("PCA top 500 variable genes")
  plot[[length(plot)+1]] <- plotParts
  
  ggsave(filename = paste(base, "/plots/QC/analysis.", round, ".PCA_noCorrection_samples.pdf", sep=""), plot = plot[[length(plot)]])
  
  #test for correlations between samples
  df2plot <- data.frame()
  for (group in unique(sample_list_xls$group)) {
    
    ids <- sample_list_xls[which(sample_list_xls$group == group), "sample_id"]
    pw_comp <- combn(ids, 2)
    
    for (i in 1: ncol(pw_comp)) {
      comp1 <- which(colnames(vsd) == pw_comp[1,i])
      comp2 <- which(colnames(vsd) == pw_comp[2,i])
      
      sample1 <- sample_list_xls[which(sample_list_xls$sample_id == pw_comp[1,i]), "sample"]
      sample2 <- sample_list_xls[which(sample_list_xls$sample_id == pw_comp[2,i]), "sample"]
      
      this_cor <- cor(vsd[,comp1], vsd[,comp2])
      
      df2plot <- rbind(df2plot, data.frame(group = group, cor=this_cor, sample1=sample1, sample2=sample2, stringsAsFactors = F))
    }
  }
  
  if (round == 1) {
    
    #PCA identifies one sample as an outlier
    idx2rem <- c(idx2rem, d[which(d$PC2 > 0 & d$age == "P0" & d$genotype == "WT"), "sample"])
    
    #sample "77860" clusters with cKO - remove
    round_rem[[2]] <- d[which(d$PC2 > 0 & d$age == "P0" & d$genotype == "WT"), "sample"] 
    
 }
  
  plot[[length(plot)+1]] <- ggplot(df2plot, aes(x=group, y=cor, group=group)) + geom_boxplot() + geom_point() + theme_classic() + ylim(0,1) + geom_hline(yintercept = 0.9)
  ggsave(filename = paste(base, "/plots/QC/correlation." , round, ".corr_pw_comp_samples.pdf", sep=""), plot = plot[[length(plot)]])
  
}

# only single sample to remove
idx2rem
# [1] "77860"

sample_list_xls <- sample_list_xls[which(!sample_list_xls$sample_id %in% idx2rem),]

# write out the final sample list with information on removed samples
sample_list_reporting <- sample_list_xls[,c("sample_id", "age", "genotype")]
sample_list_reporting$used_for_analysis <- "Y" 
sample_list_reporting <- merge(STAR_df2plot, sample_list_reporting[,c("sample_id", "used_for_analysis")], by="sample_id", all=T )
sample_list_reporting$used_for_analysis[ which( is.na(sample_list_reporting$used_for_analysis) ) ] <- "N"
sample_list_reporting <- sample_list_reporting[,c(1:4,9,11,14)]

# this file contains the combined information on sequence alignment details
wb <- createWorkbook()
addWorksheet(wb, "sample list")
writeData(wb, "sample list", sample_list_reporting)
addFilter(wb, "sample list", row = 1, cols = 1:ncol( sample_list_reporting ))
setColWidths(wb, "sample list", cols = 1:ncol(sample_list_reporting), widths="auto")
saveWorkbook(wb, file = paste(base, "Supplement/Rptor_sample_list.xlsx", sep=""), overwrite = T) 

# sample numbers:
ddply(sample_list_xls, c("age", "genotype"), summarise,
      N    = length(genotype)
)

#   age genotype N
# 1 E13      cKO 4
# 2 E13      HET 4
# 3 E13       WT 7
# 4 E16      cKO 4
# 5 E16      HET 4
# 6 E16       WT 5
# 7  P0      cKO 5
# 8  P0      HET 4
# 9  P0       WT 2

# =============================================================================
# get information on deleted region - exon6 of Rptor
# see 'del_count' folder for script on how this coverage was calculated
# files are available as Rptor_del_counts.tar.gz on GitHub
# =============================================================================

files <- list.files(path = gene_count_base, pattern="Rptor_exon6.counts.out", recursive = T)
#make path absolute
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

for(age in c("E13", "E16", "P0")) {

    #extracts samples to analyse - based an age
    sl_txt_sub <- sample_list_xls[which(sample_list_xls$age == age),]
                                        
    #extract the columns from expression data matching the list of samples extracted above
    #also extract gene_id!
    colnames_to_extract <- c("gene_id", sl_txt_sub$sample_id)
    
    #NOTE: expression cutoff used here is 10
    all_filter <- all[, which(colnames(all) %in% colnames_to_extract)]
    all_filter <- all_filter[which(rowMeans(all_filter[,2:ncol(all_filter)]) > 10),]
    
    #add delcounts data
    all_filter <- rbind(all_filter, delcounts[,colnames(all_filter)])
    
    #prepare metadata for later DESeq2 analysis    
    meta<-data.frame(  gene=sl_txt_sub$gene, genotype=sl_txt_sub$genotype, sample = sl_txt_sub$sample, age = sl_txt_sub$age, color=sl_txt_sub$color)
    rownames(meta) <- colnames(all_filter[,2:(ncol(all_filter))])
    
    #prepare countdata - change class to integer (read counts have to be integer)
    countData <- sapply(all_filter[,c(2:(ncol(all_filter)))], as.integer)
    #add rownames - name of genes
    rownames(countData) <- all_filter[,c(1)]
    
    #prepare DESeq2 comparison
    ddsMat<-DESeqDataSetFromMatrix(countData = countData,
                                   colData=meta,
                                   design= ~ genotype)
    
    #do DESeq2 analysis
    dds <- DESeq(ddsMat, fitType="local", quiet=T, betaPrior = F, parallel = T)

    #save the normalised counts
    norm_counts[[age]] <- counts(dds, normalize=T)
    
    #prepare variance stabilised counts
    vsd <- vst(dds, blind=FALSE)
    vsd_list[[age]] <- vsd
    
    #prepare a list of contrasts to be performed
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

# write out DEG analysis
# this file contains the raw DEG output
wb <- createWorkbook()
for (name in names(res)) {
  sheetName <- name
  
  addWorksheet(wb, sheetName)
  
  writeData( wb, sheetName, res[[ sheetName ]] )
  addFilter( wb, sheetName, row = 1, cols = 1:ncol(res[[ sheetName ]]) )
  setColWidths( wb, sheetName, cols = 1:ncol(res[[ sheetName ]]), widths="auto" )
}
saveWorkbook(wb, file = paste(base, "Supplement/DEG_raw_DEG_stats.xlsx", sep=""), overwrite = T) 

# Rptor exon6 (deleted exon) counts
df2plot <- data.frame()
for (age in c('E13', 'E16', 'P0')){
  tmp <- reshape2::melt(norm_counts[[age]]['Rptor_del',])
  df2plot <- rbind(df2plot, tmp)
}

df2plot$sample_id <- rownames(df2plot)
df2plot<- merge(df2plot, sample_list_xls, by='sample_id')
df2plot$genotype <- factor(df2plot$genotype, levels=c('WT', 'HET', 'cKO'))

Rptor_expr_plot <- ggplot(df2plot, aes(x=genotype, y=value)) + geom_boxplot() + geom_point() + theme_classic() + facet_grid(~age)
ggsave  ( plot = Rptor_expr_plot, filename = paste(base, "/plots/Supplement/FigS3B_Rptor_ex6_expression.pdf", sep=""), width=4, height=4)

Rptor_ex6_norm_expr <- df2plot[,c('value', 'age', 'genotype')]

# write out the raw data used for the plot
wb <- createWorkbook()
sheetName<-paste("Rptor_ex6_expression")

addWorksheet(wb, sheetName)

writeData(wb, sheetName, Rptor_ex6_norm_expr)
addFilter(wb, sheetName, row = 1, cols = 1:ncol(Rptor_ex6_norm_expr))
setColWidths(wb, sheetName, cols = 1:ncol(Rptor_ex6_norm_expr), widths="auto")

saveWorkbook(wb, file = paste(base, "Supplement/FigS3B_Rptor_ex6_expr.xlsx", sep=""), overwrite = T) 

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
#   [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggVennDiagram_1.5.2         BiocParallel_1.36.0         simplifyEnrichment_1.12.0   plyr_1.8.9                 
# [5] pheatmap_1.0.12             org.Mm.eg.db_3.18.0         AnnotationDbi_1.64.1        clusterProfiler_4.10.1     
# [9] ggplot2_3.5.0               DESeq2_1.42.1               SummarizedExperiment_1.32.0 Biobase_2.62.0             
# [13] MatrixGenerics_1.14.0       matrixStats_1.2.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.8        
# [17] IRanges_2.36.0              S4Vectors_0.40.2            BiocGenerics_0.48.1         openxlsx_4.2.5.2           
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      shape_1.4.6.1           rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3         
# [6] farver_2.1.1            ragg_1.3.0              GlobalOptions_0.1.2     fs_1.6.3                zlibbioc_1.48.2        
# [11] vctrs_0.6.5             Cairo_1.6-2             memoise_2.0.1           RCurl_1.98-1.14         ggtree_3.10.1          
# [16] S4Arrays_1.2.1          SparseArray_1.2.4       gridGraphics_0.5-1      cachem_1.0.8            igraph_2.0.2           
# [21] iterators_1.0.14        lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.6-5            R6_2.5.1               
# [26] fastmap_1.1.1           gson_0.1.0              clue_0.3-65             GenomeInfoDbData_1.2.11 digest_0.6.34          
# [31] aplot_0.2.2             enrichplot_1.22.0       colorspace_2.1-0        patchwork_1.2.0         textshaping_0.3.7      
# [36] RSQLite_2.3.6           org.Hs.eg.db_3.18.0     labeling_0.4.3          fansi_1.0.6             httr_1.4.7             
# [41] polyclip_1.10-6         abind_1.4-5             compiler_4.3.2          bit64_4.0.5             withr_3.0.0            
# [46] doParallel_1.0.17       viridis_0.6.5           DBI_1.2.2               ggforce_0.4.2           MASS_7.3-60            
# [51] proxyC_0.4.1            DelayedArray_0.28.0     rjson_0.2.23            HDO.db_0.99.1           tools_4.3.2            
# [56] ape_5.8                 scatterpie_0.2.2        zip_2.3.1               glue_1.7.0              nlme_3.1-163           
# [61] GOSemSim_2.28.1         shadowtext_0.1.3        cluster_2.1.4           reshape2_1.4.4          fgsea_1.28.0           
# [66] generics_0.1.3          gtable_0.3.4            tidyr_1.3.1             data.table_1.15.0       xml2_1.3.6             
# [71] tidygraph_1.3.1         utf8_1.2.4              XVector_0.42.0          foreach_1.5.2           ggrepel_0.9.5          
# [76] pillar_1.9.0            stringr_1.5.1           yulab.utils_0.1.4       circlize_0.4.16         splines_4.3.2          
# [81] dplyr_1.1.4             tweenr_2.0.3            treeio_1.26.0           lattice_0.21-9          bit_4.0.5              
# [86] tidyselect_1.2.1        GO.db_3.18.0            ComplexHeatmap_2.18.0   locfit_1.5-9.12         tm_0.7-14              
# [91] Biostrings_2.70.3       gridExtra_2.3           NLP_0.3-0               graphlayouts_1.1.1      stringi_1.8.3          
# [96] lazyeval_0.2.2          ggfun_0.1.4             codetools_0.2-19        ggraph_2.2.1            tibble_3.2.1           
# [101] qvalue_2.34.0           ggplotify_0.1.2         cli_3.6.2               systemfonts_1.0.6       munsell_0.5.0          
# [106] Rcpp_1.0.12             png_0.1-8               parallel_4.3.2          blob_1.2.4              DOSE_3.28.2            
# [111] bitops_1.0-7            slam_0.1-50             viridisLite_0.4.2       tidytree_0.4.6          scales_1.3.0           
# [116] packrat_0.9.2           purrr_1.0.2             crayon_1.5.2            GetoptLong_1.0.5        rlang_1.1.3            
# [121] cowplot_1.1.3           fastmatch_1.1-4         KEGGREST_1.42.0
