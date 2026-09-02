# =================================================================
# fourth step of scRNA-Seq analysis
# downstream analysis of iDEA GSEA
# Note that the actual iDEA analysis was done on the ISTA cluster
# details are in the iDEA_analysis folder
# =================================================================

library (simplifyEnrichment)
library (GOSemSim)
library (ggplot2)
library (RColorBrewer)
library (pheatmap)
library (cowplot)
library (openxlsx)

set.seed(2401)

# =======================
# set up I/O folders
# =======================
base <- "./scRNA-Seq/biorxiv/"

# lookup table from script preparing GO terms for iDEA analysis
term2id <- readRDS( file = paste(base, "/RDS_files/term2id.RDS", sep = "") )
id2term <- names( term2id )
names(id2term) <- gsub( pattern = " ", replacement = "_", x = term2id )

# prepared link between GO term and group
# NOTE: this is a file I created manually and saved
all_GO_group <- readRDS( file = paste(base, "other_data/GO_term_groups.RDS", sep=""))

# =====================
# simplifyGO analysis
# =====================

# read the iDEA results - this comes from iDEA analysis done on the ISTA cluster
# see iDEA_analysis folder for details
# Note that this file was 95M and is not provided on GitHub
idea_gen <- readRDS( paste(base, "RDS_files/E18_P1.idea_results.RDS", sep="") )

# group GO terms with with simplifyGO
iDEA_results <- list()
for (cluster in unique(names(idea_gen))) {
  iDEA_results[[ cluster ]] <- list()
  for (group in unique(names(idea_gen[[cluster]]))) {
    iDEA_results[[ cluster ]][[ group ]] <- list()
    for (dir in c("all")) {
      
      if (!any(is.na(idea_gen[[ cluster ]][[ group ]][[ dir ]]))) {
        
        tmp_df <- idea_gen[[ cluster ]][[ group ]][[ dir ]]@gsea
        
        #analyse top50 terms for each group
        tmp_df <- tmp_df[order(tmp_df$pvalue_louis, decreasing = F)[1:50],]

        terms <- id2term[ tmp_df$annot_id  ]
        
        iDEA_results[[ cluster ]][[ group ]][[ dir ]] <- tmp_df 
          
        mat = GO_similarity( terms, ont = "BP" )
        
        # simplifyGO is used to group GO terms into clusters for easier analysis
        # save the heatmap with descriptions
        pdf( file = paste(base, "Supplement/", cluster, "_", gsub(pattern = "/", replacement = "", x = group), "_", dir, "_pathway_overview_heatmap.pdf", sep="") )
        df = simplifyGO(mat, control = list("cutoff" = 0.7), order_by_size = T, verbose = F) 
        dev.off()
        
        # group the terms for easier matching of heatmap and data frame
        cluster_order <- names( sort( table(df$cluster), decreasing = T ) )
        df$ordered_cluster <- factor(df$cluster, levels = cluster_order)
        df <- df[ order(df$ordered_cluster), ]
        df$GO_term <- term2id[ df$id ]
        df$group <- all_GO_group[ gsub( pattern = " ", replacement = "_", x = df$GO_term), 2]
          
        write.table(x = df, file = paste(base, "Supplement/", cluster, "_", gsub(pattern = "/", replacement = "", x = group), "_", dir, "_pathway_overview_df.csv", sep=""), 
                       quote = F, row.names = F, sep="\t")
          
       }
    }  
  }
}

# NOTE: GO groups modified manually at this step by adding group column
#       this was in subsequent steps largely automatised by using the all_GO_group lookup table


# simple fraction of GO groups
frac_df <- data.frame()
for (name_tmp1 in c("E18_Satb2_all", "E18_Ctip2_all", "P1_Satb2_all", "P1_Ctip2_all")) {
  tmp1 <- read.table( file = paste(base, "Supplement/", name_tmp1, "_pathway_overview_df.csv", sep=""), header = T, sep = "\t")
  tmp_df <- as.data.frame( table( tmp1$group) )
  tmp_df$rel <- tmp_df$Freq / sum( tmp_df$Freq )
  tmp_df$group <- name_tmp1
  frac_df <- rbind(frac_df, tmp_df)
}

frac_df$group <- factor(frac_df$group, levels = rev(c("E18_Satb2_all", "E18_Ctip2_all", "P1_Satb2_all", "P1_Ctip2_all")))
go_group_frac <- ggplot( frac_df, aes(y=group, x=rel, fill=Var1)) + geom_bar(stat="identity", color="black") + theme_classic() + 
  scale_fill_manual( values = c("axon" = "blue", "mitochondria" = "orange", "cytokine" = "green", "apoptosis" = "grey", "ROX" = "darkviolet", "metabolism" = "white")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle("top50 GO terms")
ggsave( filename = paste(base, "plots/Main/Fig6D_go_group_frac.pdf", sep=""), plot = go_group_frac, width = 6, height = 3 ) 

# ===================================================
# prepare an overview figure for GO term enrichment
# ===================================================

# collect all GO terms in each GO term group
# these are from the top50 lists
# NOTE: there are not 50 terms in each list of the simplifyGO analysis
all_term_list <- list()
for (go_group in c("metabolism", "apoptosis", "axon", "cytokine", "mitochondria", "ROX")) {
  all_term_list[[ go_group ]] <- c()
  for (age in c("E18", "P1")) {
    for (layer in c("Ctip2", "Satb2")) {
      tmp1 <- read.table( file = paste(base, "Supplement/", age, "_", layer, "_all_pathway_overview_df.csv", sep=""), header = T, sep = "\t")
      all_term_list[[ go_group ]] <- c( all_term_list[[ go_group ]],
                                                          gsub( pattern = " ", replacement = "_", x = tmp1$GO_term[ which(tmp1$group == go_group) ] ) )
    }
  }
}

# prepare heatmaps of GO pvalue scores - one for each go group - using corrected p-values
# this is meant as a more comprehensive overview of the data

all_term_df <- list()
for (go_group in c("metabolism", "apoptosis", "axon", "cytokine", "mitochondria", "ROX")) {
  terms <- unique( all_term_list[[ go_group ]] )
  all_term_df[[ go_group ]] <- data.frame()
  
  for (age in c("E18", "P1")) {
    for (layer in c("Ctip2", "Satb2")) {
      
      tmp <- idea_gen[[ age ]][[ layer ]][[ "all" ]]@gsea
      rownames(tmp) <- tmp$annot_id
      tmp <- tmp[ terms , c( "annot_id", "pvalue_louis") ]
      colnames(tmp) <- c( "annot_id", paste(age, layer, sep="_"))
      tmp <- tmp[ which(!is.na(tmp$annot_id) ), ]
      
      if (nrow(all_term_df[[ go_group ]]) == 0){
        all_term_df[[ go_group ]] <- tmp  
      } else {
        all_term_df[[ go_group ]] <- merge( all_term_df[[ go_group ]], tmp, by="annot_id", all=T )
      }
      
    }
  }
  
  # convert the data frame to a matrix with p-value scores
  all_term_df[[ go_group ]] <- unique( all_term_df[[ go_group ]] )
  rownames( all_term_df[[ go_group ]] ) <- all_term_df[[ go_group ]]$annot_id
  all_term_df[[ go_group ]] <- all_term_df[[ go_group ]][,2:ncol(all_term_df[[ go_group ]])]
  all_term_df[[ go_group ]][ is.na(all_term_df[[ go_group ]]) ] <- 1
  all_term_df[[ go_group ]] <- log10( all_term_df[[ go_group ]] ) * -1
  all_term_df[[ go_group ]][ all_term_df[[ go_group ]] > 2 ] <- 2
  
}

# detailed analysis of specific go term groups
# for the main figure
# needs some group specific parameters for nice heat map ploting

GOgroup_parameters <- list()
GOgroup_parameters[[ "axon" ]] <- list()
GOgroup_parameters[[ "axon" ]][[1]] <- 3
GOgroup_parameters[[ "axon" ]][[2]] <- c(3,2,1)
GOgroup_parameters[[ "apoptosis" ]] <- list()
GOgroup_parameters[[ "apoptosis" ]][[1]] <- 3
GOgroup_parameters[[ "apoptosis" ]][[2]] <- c(3,1,2)
GOgroup_parameters[[ "ROX" ]] <- list()
GOgroup_parameters[[ "ROX" ]][[1]] <- 1
GOgroup_parameters[[ "ROX" ]][[2]] <- 1

term_order_list <- list()

cluster_freq_list <- list()
for (go_group in c("axon", "apoptosis", "ROX")) {
  
  # group the GO terms with an arbitrary cutoff
  clust_tree <- hclust( dist( all_term_df[[ go_group ]] ))
  clust_terms <- cutree(clust_tree, k = GOgroup_parameters[[ go_group ]][[1]] )
  
  # order them in a convenient way
  term_order_df <- data.frame( labels = clust_tree$labels[ clust_tree$order ], order =  1:length(clust_tree$order))
  
  # prepare data frames and sample order for a nice heatmap
  term_order_df$cluster <- clust_terms[ term_order_df$labels ]
  term_order_df$cluster <- factor(term_order_df$cluster, levels = GOgroup_parameters[[ go_group ]][[2]] )
  term_order_df <- term_order_df[ order(term_order_df$cluster, term_order_df$order), ]
  rownames(term_order_df) <- term_order_df$labels
  
  gaps_row <- sapply( 1:3, function (x) max(which(term_order_df$cluster == x)))
  
  pheatmap( all_term_df[[ go_group ]][ term_order_df$labels, c("E18_Satb2", "E18_Ctip2", "P1_Satb2", "P1_Ctip2")], cluster_rows = F, cluster_cols = F,
            annotation_row = term_order_df[,3,drop=F], gaps_row = gaps_row,
            filename = paste(base, "plots/Main/Fig6E_gsea.", go_group, ".score.pdf", sep="" ), width = 8)
  
  # to give an overview of the long lists of GO terms
  # calculate the frequency of each word in the GO term lists
  for (cluster in unique(term_order_df$cluster)) {
    term_parts <- sapply( term_order_df$labels[ which( term_order_df$cluster == cluster ) ], function (x) unique(strsplit(x = x, split = "_")[[1]] ) )
    term_frac <- as.data.frame( table( unlist(term_parts) ) / length( which( term_order_df$cluster == cluster ) ), decreasing = T )
    term_frac <- term_frac[ order(term_frac$Freq, decreasing = T), ]
    term_frac$Var1 <- factor( term_frac$Var1, levels = term_frac$Var1 )
    term_frac$Freq_round <- round(term_frac$Freq, 2) 
    term_frac$cluster <- cluster
    
    cluster_freq_list[[ paste(go_group, cluster, sep="_") ]] <- term_frac
    
  }
  
  # save the term order for later plotting
  term_order_list[[ go_group ]] <- term_order_df$labels
}

# save fractions of words in the GO ID
wb <- createWorkbook()

for ( name in names(cluster_freq_list)) {
  addWorksheet(wb, name)
  
  writeData(wb, name, cluster_freq_list[[ name ]])
  addFilter(wb, name, row = 1, cols = 1:ncol( cluster_freq_list[[ name ]] ))
  setColWidths(wb, name, cols = 1:ncol(cluster_freq_list[[ name ]]), widths="auto")
  
}

saveWorkbook(wb, file = paste(base, "Supplement/cluster_freq_list.xlsx", sep=""), overwrite = T) 

# prepare supplemental heatmaps to give overview of all GO groups
# and all GO terms
# letter labels of plots
letter_label <- LETTERS[1:6]
names(letter_label) <- c("apoptosis", "axon", "ROX", "cytokine", "metabolism", "mitochondria")

for (go_group in names( all_term_df )) {
  
  go_df <- all_term_df[[ go_group ]]
  
  if (is.null( term_order_list[[ go_group]] )) {
    term_order <- rownames(go_df)
    cluter_rows <- TRUE
  } else {
    term_order <- term_order_list[[ go_group]]
    cluter_rows <- FALSE
  }
  
  pheatmap(all_term_df[[ go_group ]][ term_order, c("E18_Satb2", "E18_Ctip2", "P1_Satb2", "P1_Ctip2")], cluster_rows = T, cluster_cols = F, main = go_group, 
           filename = paste(base, "plots/Supplement/FigS7", letter_label[[go_group]], "_gsea.", go_group, ".score.pdf", sep="" ), width = 10 )
  
}
         
         
# plot top 10 genes (ranked by PIP) for the gene set enrichment

for ( go_group in c("axon", "apoptosis", "ROX")) {

  # identify terms that are significantly enriched in only one comparison
  max_score <- apply(all_term_df[[ go_group ]], 1, function (x) {
    length(which(x > 1.5)) == 1 
  })
  
  # identify which comparison this term is specific for
  unique_terms <- sapply( names( which( max_score ) ), function (x) {
    which( all_term_df[[ go_group ]][x, ] == max(all_term_df[[ go_group ]][x, ]) )
  })
  
  # collect all the PIP values for these unique terms
  top_pip <- list()
  # focus only on E18
  for ( idx in c(2,1) ) {
    # extract sample name parts
    name_parts <- strsplit( colnames(all_term_df[[ go_group ]])[[idx]], "_" )[[1]]
    
    # print message for QC
    message( paste(go_group, idx, name_parts ))
             
    top_pip[[ paste(name_parts, collapse="_") ]] <- data.frame()
    
    for (id in names(which( unique_terms == idx))) {
      
      tmp_pip <- as.data.frame( idea_gen[[ name_parts[1] ]][[ name_parts[2] ]][[ "all" ]]@de[[ id ]]$pip )
      tmp_pip$gene <- rownames(tmp_pip)
      tmp_pip$id <- paste( paste(name_parts, collapse="_"), id, sep="_") 
      top_pip[[ paste(name_parts, collapse="_") ]] <- rbind( top_pip[[ paste(name_parts, collapse="_") ]], tmp_pip  )
      
    }
  }
  
  # extract the top PIP value for each gene - each GSEA analysis produced one PIP value
  top_pip_merged <- lapply( names(top_pip), function(comp_id){
    sapply( unique(top_pip[[ comp_id ]]$gene), function (gene) {
      max( top_pip[[ comp_id ]][ which( top_pip[[ comp_id ]]$gene == gene ), "PIP"] )
    })
  })
  
  # combine the top10 genes of each comparison in one figure
  pip2plot <- data.frame()
  for (i in 2:1) {
    if (length(top_pip_merged[[i]] ) > 0) {
      tmp_gl <- names( head( sort( top_pip_merged[[i]], decreasing = T), n=10) )
      
      if (length(top_pip_merged[[1]]) > 0 ) {
        tmp_df <- as.data.frame( top_pip_merged[[1]][ tmp_gl ] )
        tmp_df$gene <- tmp_gl
        tmp_df$comp <- names( top_pip )[1]
        colnames( tmp_df ) <- c( "PIP", "gene", "comp")
       } else {
        tmp_df <- data.frame( PIP = 0, gene = tmp_gl, comp = names( top_pip )[1] )
      }
      
      pip2plot <- rbind( pip2plot, tmp_df ) 
      
      if ( length (top_pip_merged[[2]]) > 0) {
        tmp_df <- as.data.frame( top_pip_merged[[2]][ tmp_gl ] )
        tmp_df$gene <- tmp_gl
        tmp_df$comp <- names( top_pip )[2]
        colnames( tmp_df ) <- c( "PIP", "gene", "comp")
      } else {
        tmp_df <- data.frame( PIP = 0, gene = tmp_gl, comp = names( top_pip )[2] )
      }
      
      pip2plot <- rbind( pip2plot, tmp_df ) 
      
    }
    
  }
  
  pip2plot$PIP[ which(is.na(pip2plot$PIP)) ] <- 0
  # create a nice order of genes to plot
  if (length(top_pip_merged[[ 2 ]]) > 0 & length(top_pip_merged[[ 1 ]]) > 0) {
    pip2plot$gene <- factor( pip2plot$gene, levels = c( names( head( sort( top_pip_merged[[ 1 ]], decreasing = T), n=10) ), 
                                                        rev ( names( head( sort( top_pip_merged[[ 2 ]], decreasing = T), n=10) ) ) ) )
    
  } else if (length(top_pip_merged[[ 1 ]]) > 0) {

    pip2plot$gene <- factor( pip2plot$gene, levels = c( names( head( sort( top_pip_merged[[ 1 ]], decreasing = T), n=10) ) ) )
    
  } else if (length(top_pip_merged[[ 2 ]]) > 0) {
    
    pip2plot$gene <- factor( pip2plot$gene, levels = c( rev( names( head( sort( top_pip_merged[[ 2 ]], decreasing = T), n=10) ) ) ) )
    
  }
  
  pip2plot$comp <- factor( pip2plot$comp, levels = c("E18_Satb2", "E18_Ctip2"))
  
  ggplot( pip2plot, aes( x=gene, y=PIP, fill=comp)) + geom_bar( stat="identity", position = "dodge", color="black" ) + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual( values = c("E18_Ctip2" = "black", "E18_Satb2" = "grey90")) + ggtitle( go_group )
  
  ggsave( filename = paste(base, "plots/Main/Fig6F_expression.", go_group, ".pdf", sep=""), width = 5, height = 3 ) 
  
}

# end of analysis

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
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] iDEA_1.0.1                openxlsx_4.2.5.2          cowplot_1.1.3             pheatmap_1.0.12           RColorBrewer_1.1-3       
# [6] ggplot2_3.5.0             GOSemSim_2.28.1           simplifyEnrichment_1.12.0 BiocGenerics_0.48.1      
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        farver_2.1.1            dplyr_1.1.4             blob_1.2.4              Biostrings_2.70.3      
# [6] bitops_1.0-7            fastmap_1.1.1           RCurl_1.98-1.14         digest_0.6.34           lifecycle_1.0.4        
# [11] cluster_2.1.4           Cairo_1.6-2             NLP_0.3-0               KEGGREST_1.42.0         RSQLite_2.3.6          
# [16] magrittr_2.0.3          compiler_4.3.2          doSNOW_1.0.20           rlang_1.1.3             tools_4.3.2            
# [21] utf8_1.2.4              labeling_0.4.3          bit_4.0.5               xml2_1.3.6              pkgload_1.3.4          
# [26] withr_3.0.0             stats4_4.3.2            fansi_1.0.6             tm_0.7-14               colorspace_2.1-0       
# [31] GO.db_3.18.0            scales_1.3.0            iterators_1.0.14        cli_3.6.2               crayon_1.5.2           
# [36] ragg_1.3.0              generics_0.1.3          rstudioapi_0.16.0       httr_1.4.7              rjson_0.2.23           
# [41] DBI_1.2.2               cachem_1.0.8            zlibbioc_1.48.2         parallel_4.3.2          AnnotationDbi_1.64.1   
# [46] XVector_0.42.0          proxyC_0.4.1            matrixStats_1.2.0       yulab.utils_0.1.4       vctrs_0.6.5            
# [51] Matrix_1.6-5            slam_0.1-50             IRanges_2.36.0          GetoptLong_1.0.5        S4Vectors_0.40.2       
# [56] bit64_4.0.5             pbmcapply_1.5.1         clue_0.3-65             systemfonts_1.0.6       packrat_0.9.2          
# [61] foreach_1.5.2           snow_0.4-4              glue_1.7.0              codetools_0.2-19        stringi_1.8.3          
# [66] shape_1.4.6.1           gtable_0.3.4            GenomeInfoDb_1.38.8     ComplexHeatmap_2.18.0   munsell_0.5.0          
# [71] tibble_3.2.1            pillar_1.9.0            GenomeInfoDbData_1.2.11 circlize_0.4.16         R6_2.5.1               
# [76] textshaping_0.3.7       doParallel_1.0.17       lattice_0.21-9          Biobase_2.62.0          png_0.1-8              
# [81] memoise_2.0.1           zip_2.3.1               Rcpp_1.0.12             org.Hs.eg.db_3.18.0     fs_1.6.3               
# [86] pkgconfig_2.0.3         GlobalOptions_0.1.2  
