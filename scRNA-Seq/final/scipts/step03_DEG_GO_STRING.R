# ==================================================================================
# use the neurons identified in an earlier step to perform DEG and network analysis
# ==================================================================================

library (Seurat)
library (ggplot2)
library (clusterProfiler)
library (dplyr)
library (tidyr)
library (openxlsx)

# group genes and display network connection of these groups
make_group_graph <- function(g, group_attr = "group", sep = ",",
                             keep_loops = FALSE) {
  
  # ------------------------------------------------------------
  # 1. Make sure the gene-level graph has vertex names
  # ------------------------------------------------------------
  # igraph edges are easiest to work with when every vertex has a
  # unique name. If the graph does not already have V(g)$name,
  # assign names using vertex indices.
  if (is.null(V(g)$name)) {
    V(g)$name <- as.character(seq_len(vcount(g)))
  }
  
  # ------------------------------------------------------------
  # 2. Define helper function to split group labels
  # ------------------------------------------------------------
  # Each gene may have a group annotation such as:
  #   "GroupA"
  # or, for multi-mapping genes:
  #   "GroupA,GroupB,GroupC"
  #
  # This function splits such strings into a character vector:
  #   c("GroupA", "GroupB", "GroupC")
  #
  # It also removes leading/trailing spaces and removes duplicated
  # group names within the same gene.
  split_groups <- function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return(character(0))
    unique(trimws(strsplit(as.character(x), sep, fixed = TRUE)[[1]]))
  }
  
  # ------------------------------------------------------------
  # 3. Extract and split the group annotation for every gene/node
  # ------------------------------------------------------------
  # vertex_attr(g, group_attr) extracts the vertex attribute containing
  # the group labels. For example, if group_attr = "group", this extracts
  # V(g)$group.
  #
  # groups_list is a list where each element contains the group or groups
  # assigned to one gene.
  groups_list <- lapply(igraph::vertex_attr(g, group_attr), split_groups)
  
  # ------------------------------------------------------------
  # 4. Create a gene-level table
  # ------------------------------------------------------------
  # This table has one row per gene.
  #
  # gene   = gene/node name
  # groups = list-column containing all group assignments for that gene
  # nmap   = number of groups that gene maps to
  #
  # Genes with no group assignment are removed.
  gene_tbl <- tibble(
    gene = V(g)$name,
    groups = groups_list,
    nmap = lengths(groups_list)
  ) %>%
    filter(nmap > 0)
  
  # ------------------------------------------------------------
  # 5. Expand gene-to-group mapping
  # ------------------------------------------------------------
  # Convert the list-column into a long table.
  #
  # Example:
  #   gene1 -> c("A", "B")
  #
  # becomes:
  #   gene1   A
  #   gene1   B
  #
  # This is needed because multi-mapping genes must contribute to
  # multiple group-level nodes.
  gene_group <- gene_tbl %>%
    unnest_longer(groups, values_to = "group")
  
  # ------------------------------------------------------------
  # 6. Create the group-level vertex table
  # ------------------------------------------------------------
  # The resulting group network will have one node per group.
  #
  # n_genes:
  #   Counts each gene once in every group it belongs to.
  #   For example, if gene1 belongs to groups A and B, it adds 1
  #   to both A and B.
  #
  # n_genes_fractional:
  #   Fractional count adjusted for multi-mapping.
  #   For example, if gene1 belongs to two groups, it contributes
  #   1/2 to each group. If it belongs to three groups, it contributes
  #   1/3 to each group.
  #
  # The final vertex table must have a column named "name" because
  # graph_from_data_frame() uses that column to assign vertex names.
  
  total_genes <- igraph::vcount(g)
  
  group_vertices <- gene_group %>%
    group_by(group) %>%
    summarise(
      n_genes = n_distinct(gene),
      n_genes_fractional = sum(1 / nmap),
      .groups = "drop"
    ) %>%
    mutate(
      n_genes_fractional_percent = 100 * n_genes_fractional / total_genes
    ) %>%
    rename(name = group)
  
  # ------------------------------------------------------------
  # 7. Extract the original gene-gene edge list
  # ------------------------------------------------------------
  # ends(g, E(g), names = TRUE) returns a two-column matrix containing
  # the names of the two genes connected by each edge.
  #
  # edge_id is added so we can track each original edge after expansion
  # to group-group combinations.
  ee <- igraph::ends(g, igraph::E(g), names = TRUE)
  
  edge_tbl <- tibble(
    edge_id = seq_len(igraph::ecount(g)),
    gene_from = ee[, 1],
    gene_to   = ee[, 2]
  )
  
  # ------------------------------------------------------------
  # 8. Expand gene-gene edges to group-group edges
  # ------------------------------------------------------------
  # Each gene-gene edge is replaced by all possible combinations of
  # groups assigned to the two genes.
  #
  # Example:
  #   gene1 belongs to A,B
  #   gene2 belongs to C
  #   gene1 -- gene2
  #
  # becomes:
  #   A -- C
  #   B -- C
  #
  # If both genes map to multiple groups, all combinations are created.
  #
  # nmap_from and nmap_to are retained so fractional edge weights can
  # be calculated later.
  edge_groups <- edge_tbl %>%
    left_join(
      gene_group %>%
        select(gene_from = gene, group_from = group, nmap_from = nmap),
      by = "gene_from"
    ) %>%
    left_join(
      gene_group %>%
        select(gene_to = gene, group_to = group, nmap_to = nmap),
      by = "gene_to"
    ) %>%
    filter(!is.na(group_from), !is.na(group_to))
  
  # ------------------------------------------------------------
  # 9. Define group-level edge direction
  # ------------------------------------------------------------
  # For directed graphs:
  #   Preserve the direction from group_from to group_to.
  #
  # For undirected graphs:
  #   Sort the group names so that A--B and B--A are treated as the
  #   same edge.
  if (igraph::is_directed(g)) {
    edge_groups <- edge_groups %>%
      mutate(
        from = group_from,
        to   = group_to
      )
  } else {
    edge_groups <- edge_groups %>%
      mutate(
        from = pmin(group_from, group_to),
        to   = pmax(group_from, group_to)
      )
  }
  
  # ------------------------------------------------------------
  # 10. Optionally remove within-group edges
  # ------------------------------------------------------------
  # A loop is an edge where from == to.
  #
  # In this context, loops represent original gene-gene edges where
  # both genes share the same group.
  #
  # If keep_loops = FALSE, those within-group connections are removed.
  # If keep_loops = TRUE, they are retained as self-edges in the
  # group-level graph.
  if (!keep_loops) {
    edge_groups <- edge_groups %>%
      filter(from != to)
  }
  
  # ------------------------------------------------------------
  # 11. Summarise expanded edges into group-group edges
  # ------------------------------------------------------------
  # After expansion, there may be many rows for the same group pair.
  # These are collapsed into a single edge with weights.
  #
  # weight:
  #   Number of distinct original gene-gene edges connecting the two
  #   groups. This prevents a single original edge from being counted
  #   multiple times for the same group pair.
  #
  # weight_expanded:
  #   Number of expanded group-group connections. This counts every
  #   group-combination generated from multi-mapping genes.
  #
  # weight_fractional:
  #   Fractional edge count adjusted for multi-mapping genes.
  #   If one gene belongs to 2 groups and the other belongs to 3 groups,
  #   that gene-gene edge contributes 1 / (2 * 3) to each expanded
  #   group-group edge.
  # Total number of edges in the original gene-level network
  total_edges <- igraph::ecount(g)
  
  group_edges <- edge_groups %>%
    group_by(from, to) %>%
    summarise(
      # Number of distinct original gene-gene edges connecting these groups
      weight = n_distinct(edge_id),
      
      # Number of expanded group-group edges after accounting for multi-mapping
      weight_expanded = n(),
      
      # Fractional edge weight adjusted for multi-mapping genes
      weight_fractional = sum(1 / (nmap_from * nmap_to)),
      
      .groups = "drop"
    ) %>%
    mutate(
      # Raw edge count relative to total number of original edges
      weight_relative = weight / total_edges,
      
      # Fractional edge count relative to total number of original edges
      weight_fractional_percent = 100 * weight_fractional / total_edges,
      
    )
  
  # ------------------------------------------------------------
  # 12. Build the final group-level igraph object
  # ------------------------------------------------------------
  # graph_from_data_frame() creates an igraph object from:
  #
  # d:
  #   Edge table. Must contain columns called "from" and "to".
  #
  # vertices:
  #   Vertex table. Must contain a column called "name".
  #
  # directed:
  #   Uses the same directed/undirected status as the original graph.
  group_graph <- igraph::graph_from_data_frame(
    d = group_edges,
    vertices = group_vertices,
    directed = igraph::is_directed(g)
  )
  
  # ------------------------------------------------------------
  # 13. Return the final igraph object
  # ------------------------------------------------------------
  return(group_graph)
}


set.seed(2401)

# =======================
# set up I/O folders
# =======================
old_base <- "./scRNA-Seq/biorxiv/"
base <- "./scRNA-Seq/final/"
bulk_base <- "./bulk_RNASeq/final/" # main working folder bulk data

# ==============
# DEG analysis
# ==============

# read the neurons identified earlier
neuron_seurat <- readRDS( file = paste(base, "/RDS_files/neuron_seurat.RDS", sep=""))

# create a genotype variable that is common to all ages
neuron_seurat@meta.data$genotype <- sapply( neuron_seurat@meta.data$genotype, function(x) strsplit(x, "_")[[1]][2] )

# DEG analysis
Idents( neuron_seurat ) <- "genotype"
neuron_seurat <- JoinLayers( neuron_seurat )
deg_analysis <- FindMarkers( neuron_seurat, ident.1 = "Het", ident.2 = "WT", test.use = "MAST")
deg_analysis$symbol <- rownames(deg_analysis)

# simple plot of DEG number
deg2plot <- data.frame( n = c( length( which(deg_analysis$p_val_adj < 0.05 & deg_analysis$avg_log2FC > 0)),
length( which(deg_analysis$p_val_adj < 0.05 & deg_analysis$avg_log2FC < 0) ) ), dir = c("up", "down") )

#  n  dir
# 33   up
# 75 down

ggplot( deg2plot, aes(x=1, y=n, fill=dir) ) + geom_bar(stat = "identity") + theme_classic() +
  scale_fill_manual( values = c("up" = "gold", "down" = "cornflowerblue")) + ylab("# DEGs") + xlab("")
ggsave( paste(base, "plots/Supplement/FigS9b_DEG_number_bar.pdf", sep=""), width = 2)

# read GO terms
# this is the same as for the bulk analysis
GO_GENE_DF <- readRDS(  paste(bulk_base, "/RDS_files/GO_Term_df_collapsed.RDS", sep="")  )
colnames( GO_GENE_DF ) <- c("TERM", "GENE")
# give some group a simpler name
GO_GENE_DF$TERM[ which(GO_GENE_DF$TERM == "reactive oxygen")] <- "ROS"
GO_GENE_DF$TERM[ which(GO_GENE_DF$TERM == "epidermal growth factor")] <- "Egfr"
GO_GENE_DF$TERM[ which(GO_GENE_DF$TERM == "mitochondr")] <- "mito"

GO_GENE_DF <- GO_GENE_DF[ which(!GO_GENE_DF$TERM %in% c("TOR/PI3K", "MAPK", "Egfr", "Notch", "JAK-STAT", "mitotic cell cycle")), ] 

# ============
# GO analysis
# ============

eGO <- enricher ( gene = rownames(deg_analysis)[which(deg_analysis$p_val_adj < 0.05)], 
                      TERM2GENE = GO_GENE_DF[,c("TERM", "GENE")], 
                      pvalueCutoff = 0.99, qvalueCutoff = 0.99, 
                      pAdjustMethod = "BH",
                      minGSSize = 10, 
                      maxGSSize = 3000 )

# write out DEG and GO data
wb <- createWorkbook()

sheetName <- "DEG Het vs Ctrl"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, deg_analysis )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( deg_analysis ) )
setColWidths( wb, sheetName, cols = 1:ncol( deg_analysis ), widths="auto" )

sheetName <- "GO analysis"
addWorksheet(wb, sheetName)

writeData( wb, sheetName, eGO@result )
addFilter( wb, sheetName, row = 1, cols = 1:ncol( eGO@result ) )
setColWidths( wb, sheetName, cols = 1:ncol( eGO@result ), widths="auto" )

saveWorkbook(wb, file = paste(base, "/Supplement/DEG_GO.xlsx", sep=""), overwrite = T) 

# combine the genes into groups based on GO analysis
combined_genes <- eGO@result %>%
  select(ID, geneID) %>%
  bind_rows() %>%
  separate_rows(geneID, sep = "/") %>%
  rename(gene = geneID, Description = ID) %>%
  distinct(gene, Description) %>%
  group_by(gene) %>%
  summarise(
    Description = paste(unique(Description), collapse = ", "),
    .groups = "drop"
  )

# write out this data to use for STRING analysis
write.csv( x = combined_genes, file = paste(base, "/Supplement/DEGs_4nw.csv", sep="") )

# =====================
# network analysis
# prerequisite is the STRING network analysis on https://string-db.org/ v12.0
# =====================

  # analyse STRINGdb networks - files were downloaded from https://string-db.org/ analysis
  # setting of this analysis are documented here: Supplement/STRING_settings.png
  # stats of this network are documented here: Supplement/STRING_analysis.png
  # mapping of genes to string is documented here: Supplement/DEG_string_mapping.tsv
  # gene-gene interactions from database: string_interactions_short.tsv

library (igraph)
library (RCy3)

# analyse STRINGdb networks
complete_network <- read.table( file = paste(base, "/Supplement/string_interactions_short.tsv", sep=""))

g <- graph_from_data_frame(complete_network[,1:2], directed=F)

#remove multiple edges
g_simple <- delete_edges(g, E(g)[which_multiple(g)])

g_simple_clusters <- components(g_simple)

#get the largest connected cluster in the network
idx <- which(g_simple_clusters$csize == max(g_simple_clusters$csize))
connected_genes <- names(which(g_simple_clusters$membership == idx))

largest_connected_nw <- complete_network[ which( complete_network$V1 %in% connected_genes | complete_network$V2 %in% connected_genes ), ]
g <- graph_from_data_frame(largest_connected_nw, directed=F)

#remove multiple edges
g_simple <- delete_edges(g, E(g)[which_multiple(g)])

# read meta data from GO analysis
gene2go <- combined_genes$Description
names(gene2go) <- combined_genes$gene

# add group variable to network
V(g_simple)$group <- gene2go[ V(g_simple)$name ]

# record network sizes
length(g_simple)
# 33

# create a simplified network only showing functional groups of genes
g_group <- make_group_graph(g_simple, group_attr = "group")

# send this network to cytoscape for visualisation
createNetworkFromIgraph(g_group, title = "single_cell_v1", collection = "grouped_networks")

# ====================================================================================
# the resulting Cytoscape network was saves as Supplement/group_network.cys
# the figure was exported from Cytoscape: plots/Supplement/FigS9c_single_cell_v1.pdf
# ====================================================================================

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
#   [1] RCy3_2.22.1            igraph_2.1.4           openxlsx_4.2.5.2       tidyr_1.3.1            dplyr_1.1.4           
# [6] clusterProfiler_4.10.1 ggplot2_3.5.2          Seurat_5.0.1           SeuratObject_5.0.1     sp_2.1-3              
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                    matrixStats_1.2.0           spatstat.sparse_3.0-3       bitops_1.0-7               
# [5] enrichplot_1.22.0           HDO.db_0.99.1               httr_1.4.7                  RColorBrewer_1.1-3         
# [9] repr_1.1.7                  tools_4.3.2                 sctransform_0.4.1           backports_1.4.1            
# [13] utf8_1.2.4                  R6_2.5.1                    lazyeval_0.2.2              uwot_0.1.16                
# [17] withr_3.0.0                 prettyunits_1.2.0           gridExtra_2.3               progressr_0.14.0           
# [21] cli_3.6.2                   Biobase_2.62.0              spatstat.explore_3.2-6      fastDummies_1.7.3          
# [25] scatterpie_0.2.2            spatstat.data_3.0-4         ggridges_0.5.6              pbapply_1.7-2              
# [29] yulab.utils_0.1.4           pbdZMQ_0.3-14               gson_0.1.0                  DOSE_3.28.2                
# [33] parallelly_1.37.0           rstudioapi_0.16.0           RSQLite_2.3.6               generics_0.1.3             
# [37] gridGraphics_0.5-1          gtools_3.9.5                ica_1.0-3                   spatstat.random_3.2-2      
# [41] zip_2.3.1                   GO.db_3.18.0                Matrix_1.6-5                fansi_1.0.6                
# [45] S4Vectors_0.40.2            abind_1.4-5                 lifecycle_1.0.4             SummarizedExperiment_1.32.0
# [49] gplots_3.1.3.1              qvalue_2.34.0               SparseArray_1.2.4           Rtsne_0.17                 
# [53] grid_4.3.2                  blob_1.2.4                  promises_1.2.1              crayon_1.5.2               
# [57] miniUI_0.1.1.1              lattice_0.21-9              cowplot_1.1.3               KEGGREST_1.42.0            
# [61] pillar_1.9.0                fgsea_1.28.0                GenomicRanges_1.54.1        future.apply_1.11.1        
# [65] codetools_0.2-19            fastmatch_1.1-4             leiden_0.4.3.1              glue_1.7.0                 
# [69] packrat_0.9.2               ggfun_0.1.4                 data.table_1.15.0           vctrs_0.6.5                
# [73] png_0.1-8                   treeio_1.26.0               spam_2.10-0                 gtable_0.3.6               
# [77] cachem_1.0.8                S4Arrays_1.2.1              mime_0.12                   tidygraph_1.3.1            
# [81] survival_3.5-7              SingleCellExperiment_1.24.0 ellipsis_0.3.2              fitdistrplus_1.1-11        
# [85] ROCR_1.0-11                 nlme_3.1-163                ggtree_3.10.1               bit64_4.0.5                
# [89] progress_1.2.3              RcppAnnoy_0.0.22            GenomeInfoDb_1.38.8         irlba_2.3.5.1              
# [93] KernSmooth_2.23-22          colorspace_2.1-0            BiocGenerics_0.48.1         DBI_1.2.2                  
# [97] tidyselect_1.2.1            curl_5.2.0                  bit_4.0.5                   compiler_4.3.2             
# [101] graph_1.80.0                DelayedArray_0.28.0         plotly_4.10.4               shadowtext_0.1.3           
# [105] scales_1.4.0                caTools_1.18.2              lmtest_0.9-40               stringr_1.5.1              
# [109] digest_0.6.34               goftest_1.2-3               spatstat.utils_3.1-0        XVector_0.42.0             
# [113] htmltools_0.5.7             pkgconfig_2.0.3             base64enc_0.1-3             MatrixGenerics_1.14.0      
# [117] fastmap_1.1.1               rlang_1.1.3                 htmlwidgets_1.6.4           shiny_1.8.0                
# [121] farver_2.1.1                zoo_1.8-12                  jsonlite_1.8.8              BiocParallel_1.36.0        
# [125] GOSemSim_2.28.1             RCurl_1.98-1.14             magrittr_2.0.3              GenomeInfoDbData_1.2.11    
# [129] ggplotify_0.1.2             dotCall64_1.1-1             patchwork_1.2.0             IRkernel_1.3.2             
# [133] Rcpp_1.0.12                 ape_5.8                     viridis_0.6.5               reticulate_1.44.1          
# [137] stringi_1.8.3               ggraph_2.2.1                RJSONIO_2.0.0               zlibbioc_1.48.2            
# [141] MASS_7.3-60                 plyr_1.8.9                  MAST_1.28.0                 parallel_4.3.2             
# [145] listenv_0.9.1               ggrepel_0.9.5               deldir_2.0-2                Biostrings_2.70.3          
# [149] graphlayouts_1.1.1          IRdisplay_1.1               splines_4.3.2               tensor_1.5                 
# [153] hms_1.1.3                   uuid_1.2-0                  spatstat.geom_3.2-8         base64url_1.4              
# [157] RcppHNSW_0.6.0              reshape2_1.4.4              stats4_4.3.2                XML_3.99-0.16.1            
# [161] evaluate_0.23               tweenr_2.0.3                httpuv_1.6.14               RANN_2.6.1                 
# [165] purrr_1.0.2                 polyclip_1.10-6             future_1.33.1               scattermore_1.2            
# [169] ggforce_0.4.2               xtable_1.8-4                RSpectra_0.16-1             tidytree_0.4.6             
# [173] later_1.3.2                 viridisLite_0.4.2           tibble_3.2.1                aplot_0.2.2                
# [177] memoise_2.0.1               AnnotationDbi_1.64.1        IRanges_2.36.0              cluster_2.1.4              
# [181] globals_0.16.2
