# =======================================================================================
# step 3: of revision re-analysis
# display genes as STRING networks
# prerequisite is the STRING network analysis on https://string-db.org/ v12.0
# =======================================================================================

library (igraph)
library (RCy3)
library (dplyr)
library (tidyr)

# ==================
# helper functions 
# ==================

# copy the node layout from one network onto another
copyNodePositions <- function(source.net, target.net) {
  src.nodes <- getAllNodes(network = source.net)
  tgt.nodes <- getAllNodes(network = target.net)
  
  pos <- getNodePosition(src.nodes, network = source.net)
  pos$name <- rownames(pos)
  
  # Try to identify column names
  name.col <- intersect(c("node.names", "name", "node"), colnames(pos))[1]
  x.col    <- intersect(c("x_location", "x", "NODE_X_LOCATION"), colnames(pos))[1]
  y.col    <- intersect(c("y_location", "y", "NODE_Y_LOCATION"), colnames(pos))[1]
  
  if (is.na(name.col) || is.na(x.col) || is.na(y.col)) {
    stop("Could not identify node name/x/y columns. Check colnames(pos).")
  }
  
  pos.common <- pos[pos[[name.col]] %in% tgt.nodes, ]
  
  setNodePositionBypass(
    node.names  = pos.common[[name.col]],
    new.x.locations = pos.common[[x.col]],
    new.y.locations = pos.common[[y.col]],
    network     = target.net
  )
}


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

# check whether RCy3 could connect to cytoscape
cytoscapePing()

# =======================
# set up I/O folders
# =======================
bulk_base <- "./bulk_RNASeq/biorxiv/" # initial analysis
main_base <- "./bulk_RNASeq/final/" # main working folder

# container for network stats
network_stats <- data.frame()

for ( prefix in c( "cKO", "het" )) {
  
  # analyse STRINGdb networks - files were downloaded from https://string-db.org/ analysis
  # setting of this analysis are documented here: other_data/[prefix]_stringDB_settings.png
  # stats of this network are documented here: other_data/[prefix]_stringDB_analysis.png
  # mapping of genes to string is documented here: other_data/[prefix]_string_mapping.tsv
  
  complete_network <- read.table( file = paste(main_base, "/other_data/", prefix, "_string_interactions_short.tsv", sep=""))
  
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
  combined_genes <- read.csv( file = paste(main_base, "/Supplement/", prefix, "_genes_4nw.csv", sep=""), row.names = 1 )
  gene2go <- combined_genes$Description
  names(gene2go) <- combined_genes$gene
  
  # add group variable to network
  V(g_simple)$group <- gene2go[ V(g_simple)$name ]
  
  # record network sizes
  network_stats <- rbind( network_stats, data.frame( comparison = prefix, size = length(g_simple) ) )
  
  # create a simplified network only showing functional groups of genes
  g_group <- make_group_graph(g_simple, group_attr = "group")
  
  # send this network to cytoscape for visualisation
  createNetworkFromIgraph(g_group, title = paste(prefix, "network_v1"), collection = "grouped_networks")
}

network_stats
#   comparison size
# 1        cKO  433
# 2        het  210

# optional
# copy layout from one network to the other to allow direct comparison
copyNodePositions("het network_v1", "cKO network_v1")


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
#   [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=de_AT.UTF-8       
# [4] LC_COLLATE=en_GB.UTF-8     LC_MONETARY=de_AT.UTF-8    LC_MESSAGES=en_GB.UTF-8   
# [7] LC_PAPER=de_AT.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C             LC_MEASUREMENT=de_AT.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Vienna
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] tidyr_1.3.1  dplyr_1.1.4  RCy3_2.22.1  igraph_2.1.4
# 
# loaded via a namespace (and not attached):
#   [1] base64url_1.4       jsonlite_1.8.8      compiler_4.3.2      crayon_1.5.2       
# [5] gtools_3.9.5        tidyselect_1.2.1    IRdisplay_1.1       bitops_1.0-7       
# [9] uuid_1.2-0          fastmap_1.1.1       IRkernel_1.3.2      R6_2.5.1           
# [13] generics_0.1.3      curl_5.2.0          BiocGenerics_0.48.1 backports_1.4.1    
# [17] XML_3.99-0.16.1     graph_1.80.0        tibble_3.2.1        pillar_1.9.0       
# [21] RColorBrewer_1.1-3  rlang_1.1.3         utf8_1.2.4          stringi_1.8.3      
# [25] repr_1.1.7          fs_1.6.3            caTools_1.18.2      RJSONIO_2.0.0      
# [29] cli_3.6.2           withr_3.0.0         magrittr_2.0.3      digest_0.6.34      
# [33] rstudioapi_0.16.0   base64enc_0.1-3     packrat_0.9.2       pbdZMQ_0.3-14      
# [37] lifecycle_1.0.4     vctrs_0.6.5         KernSmooth_2.23-22  evaluate_0.23      
# [41] glue_1.7.0          stats4_4.3.2        RCurl_1.98-1.14     fansi_1.0.6        
# [45] purrr_1.0.2         httr_1.4.7          tools_4.3.2         pkgconfig_2.0.3    
# [49] htmltools_0.5.7     gplots_3.1.3.1    
