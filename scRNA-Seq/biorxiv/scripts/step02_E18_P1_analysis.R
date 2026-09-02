# ===================================
# second step of the scRNA-Seq analysis
# main analysis of E18/P1 data
# ===================================

library (clusterProfiler)
library (org.Mm.eg.db)
library (Seurat)
library (plyr)
library (ggplot2)
library (GeneOverlap)
library (GenBinomApps)

set.seed(2401)

# =======================
# set up I/O folders
# =======================
base <- "./scRNA-Seq/biorxiv/"
noSaveBase <- "./Villalba_et_al_noSave" # folder for large files that are not necessary to backup

# ======================================
# prepare reference for label transfer
# ======================================
seurat.obj.lst <- readRDS( file = paste( base, "RDS_files/DiBella.seurat.obj.list.RDS", sep="") )

# NOTE: I only use 1 P1 data set as reference
reference <- seurat.obj.lst[[1]]
# standard workflow and some QC
reference <- NormalizeData(reference, normalization.method = "LogNormalize", scale.factor = 10000)
reference <- FindVariableFeatures( reference, selection.method = "vst", nfeatures = 2000 )
reference <- ScaleData( reference )
reference <- RunPCA( reference )

ElbowPlot( reference )
# 15 seems a reasonable number

reference <- RunUMAP( reference, dims = 1:15 )
DimPlot( reference, reduction = "umap", group.by = "orig.ident")
ggsave( paste(base, "plots/QC/Reference_P1_umap.pdf", sep=""))
DimPlot( reference, reduction = "umap", group.by = "New_cellType", label=T)
ggsave( paste(base, "plots/QC/Reference_P1_umap_cellTypes.pdf", sep=""), width=10, height = 7)

# prepare cell cycle genes
cell_cycle_phase_gene <- read.table(file = paste(base, "other_data/Mus_musculus.csv", sep=""), sep = ",", header = T)
gene_conversion <- bitr(geneID = cell_cycle_phase_gene$geneID, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)
cell_cycle_phase_gene <- merge(cell_cycle_phase_gene, gene_conversion, by.x="geneID", by.y="ENSEMBL")

# ======================
# prepare our own data
# ======================

# Note that filtered_feature_bc_matrix.h5 as well as HTO_demux.csv files are available on GEO

count_matrix_E18 <- Read10X_h5( filename = paste(noSaveBase, "E18/filtered_feature_bc_matrix.h5", sep="") )
count_matrix_P1 <- Read10X_h5( filename = paste(noSaveBase, "P1/filtered_feature_bc_matrix.h5", sep="") )

seurat.object <- list()
seurat.object[["E18"]] <- CreateSeuratObject(counts = count_matrix_E18[[1]], min.cells = 3, min.features = 500)
seurat.object[["P1"]] <- CreateSeuratObject(counts = count_matrix_P1[[1]], min.cells = 3, min.features = 500)

# add genotype information
# Note that Rptor-MADM and Control-MADM samples were labeled with CMOs and annalysed together
# CMO assignment here

cmo2genotype <- c()
cmo2genotype <- c("P1_WT", "P1_Het", "E18_WT", "E18_Het", "Doublet", "Negative")
names( cmo2genotype ) <- c("HTO-CMO_1", "HTO-CMO_2", "HTO-CMO_3", "HTO-CMO_4", "Doublet", "Negative")

# add data on CMO assignment
meta_data <- read.csv( file = paste(noSaveBase, "E18/HTO_demux.csv", sep="") )
rownames( meta_data ) <- paste(meta_data$index, "1", sep="-")
seurat.object[["E18"]] <- AddMetaData( object = seurat.object[["E18"]], metadata = meta_data )
seurat.object[["E18"]]$age <- "E18" 
seurat.object[["E18"]]@meta.data$genotype <- cmo2genotype[ seurat.object[["E18"]]$hto_demux ]
seurat.object[["E18"]]@meta.data$genotype <- factor( seurat.object[["E18"]]@meta.data$genotype, 
                                                     levels = c("E18_WT", "E18_Het", "Doublet", "Negative"))

meta_data <- read.csv( file = paste(noSaveBase, "P1/HTO_demux.csv", sep="") )
rownames( meta_data ) <- paste(meta_data$index, "1", sep="-")
seurat.object[["P1"]] <- AddMetaData( object = seurat.object[["P1"]], metadata = meta_data )
seurat.object[["P1"]]$age <- "P1"
seurat.object[["P1"]]@meta.data$genotype <- cmo2genotype[ seurat.object[["P1"]]$hto_demux ]
seurat.object[["P1"]]@meta.data$genotype <- factor( seurat.object[["P1"]]@meta.data$genotype,
                                                    levels = c("P1_WT", "P1_Het", "Doublet", "Negative"))

# filter, correct for cell cycle and do label transfer
for (i in 1:length(seurat.object)) {
  
  seurat.object[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat.object[[i]], pattern = "^mt-")
  seurat.object[[i]] <- subset(seurat.object[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & nCount_RNA < 20000 & percent.mt < 5)
  seurat.object[[i]] <- NormalizeData(seurat.object[[i]], verbose = FALSE)
  seurat.object[[i]] <- FindVariableFeatures(seurat.object[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  seurat.object[[i]] <- CellCycleScoring(
    object = seurat.object[[i]],
    g2m.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "G2/M"), "SYMBOL"],
    s.features = cell_cycle_phase_gene[which(cell_cycle_phase_gene$phase == "S"), "SYMBOL"]
  )
  
  # partially regress out cell cycle effects
  seurat.object[[i]]$CC.Difference <- seurat.object[[i]]$S.Score - seurat.object[[i]]$G2M.Score
  seurat.object[[i]] <- ScaleData( seurat.object[[i]], vars.to.regress = "CC.Difference" )
  
  seurat.object[[i]] <- RunPCA(seurat.object[[i]], npcs = 30, verbose = FALSE)
  seurat.object[[i]] <- RunUMAP(seurat.object[[i]], reduction = "pca", dims = 1:30)
  seurat.object[[i]] <- FindNeighbors(seurat.object[[i]], reduction = "pca", dims = 1:30)
  seurat.object[[i]] <- FindClusters(seurat.object[[i]], resolution = 1)
  
  # Label transfer from diBella
  anchors <- FindTransferAnchors(reference = reference, query = seurat.object[[i]], dims = 1:20,
                                 reference.reduction = "pca" )
  
  predictions <- TransferData(anchorset = anchors, refdata = reference$New_cellType, dims = 1:15)
  seurat.object[[i]] <- AddMetaData(seurat.object[[i]], metadata = predictions)
  
}

# basic stats for each object
seurat.object[[1]]

# An object of class Seurat 
# 20419 features across 10612 samples within 1 assay 
# Active assay: RNA (20419 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

seurat.object[[2]]

# An object of class Seurat 
# 20851 features across 9365 samples within 1 assay 
# Active assay: RNA (20851 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

# =============
##  E18
# =============

# some plots to visualize the data
plot1 <- DimPlot( seurat.object[[ 1 ]], reduction = "umap", group.by = "predicted.id", label.size = 7, label = F) + ggtitle("E18")
plot2 <- DimPlot( seurat.object[[ 1 ]], reduction = "umap", group.by = "seurat_clusters", label.size = 7,label = T) + NoLegend()

plot1 + plot2
ggsave( paste(base, "/plots/Supplement/FigS6BC_E18_cell_type_assignment_clustering.pdf", sep=""), width=14, height = 7)

# main figure - color the neurons
seurat.object[[ 1 ]]@meta.data$main_figure_coloring <- "other"
seurat.object[[ 1 ]]@meta.data$main_figure_coloring[ which( seurat.object[[ 1 ]]@meta.data$seurat_clusters %in% c("1", "9", "7", "3", "2", "0" ))] <- "neuron"

# main figure is cropped using ylim!
DimPlot( seurat.object[[ 1 ]], reduction = "umap", group.by = "main_figure_coloring", label.size = 7,label = F) + 
  ylim(-10,5) + scale_color_manual( values = c("other" = "grey80", "neuron" = "black")) +
  ggtitle( "E18") 

# Warning message:
#   Removed 179 rows containing missing values or values outside the scale range (`geom_point()`). 

ggsave( paste(base, "/plots/Main/Fig6B_E18_all_cells_UMAP.pdf", sep=""), width=5, height = 5)

# marker expression on UMAP
marker_plot <- FeaturePlot( seurat.object[[ 1 ]], features = c("Hes5", "Mki67", "Dcx", "Eomes", "Neurod1", "Gad2", "Aldh1l1", "Sox9", "Sox10"), order = T, min.cutoff = "q15" ) & NoLegend()

marker_plot <- marker_plot & NoLegend() + theme( axis.line=element_blank(),axis.text.x=element_blank(),
                      axis.text.y=element_blank(),axis.ticks=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank() )

ggsave( filename = paste(base, "/plots/Supplement/FigS6D_E18_marker_expression_plots_noaxis.pdf", sep=""), plot = marker_plot, width=9, height=9 )

# basic QC of data set
VlnPlot(object = seurat.object[[ 1 ]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "genotype")
ggsave( paste(base, "/plots/Supplement/FigS6A_E18_basic_QC.pdf", sep=""), width=7, height=4)

# =======
##  P1
# =======
plot1 <- DimPlot( seurat.object[[ 2 ]], reduction = "umap", group.by = "predicted.id", label.size = 7, label = F) + ggtitle("P1")
plot2 <- DimPlot( seurat.object[[ 2 ]], reduction = "umap", group.by = "seurat_clusters", label.size = 7,label = T) + NoLegend()
plot1 + plot2
ggsave( paste(base, "/plots/Supplement/FigS6FG_P1_cell_type_assignment_clustering.pdf", sep=""), width=14, height = 7)

# for main figure show part of the UMAP
seurat.object[[ 2 ]]@meta.data$main_figure_coloring <- "other"
seurat.object[[ 2 ]]@meta.data$main_figure_coloring[ which( seurat.object[[ 2 ]]@meta.data$seurat_clusters %in% c(4,0))] <- "neuron"

# main figure is cropped using ylim!

DimPlot( seurat.object[[ 2 ]], reduction = "umap", group.by = "main_figure_coloring", label.size = 7,label = F) + 
  ylim(-10,5) + scale_color_manual( values = c("other" = "grey80", "neuron" = "black")) +
  ggtitle( "P1" ) + NoLegend()

# Warning message:
#   Removed 389 rows containing missing values or values outside the scale range (`geom_point()`). 

ggsave( paste(base, "/plots/Main/Fig6B_P1_all_cells_UMAP.pdf", sep=""), width=5, height = 5)

# marker expression on UMAP
marker_plot <- FeaturePlot( seurat.object[[ 2 ]], features = c("Hes5", "Mki67", "Dcx", "Eomes", "Neurod1", "Gad2", "Aldh1l1", "Sox9", "Sox10"), order = T, min.cutoff = "q15" ) & NoLegend()

marker_plot <- marker_plot & NoLegend() + theme( axis.line=element_blank(),axis.text.x=element_blank(),
                                                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                                                 axis.title.x=element_blank(),
                                                 axis.title.y=element_blank() )

ggsave( filename = paste(base, "/plots/Supplement/FigS6H_P1_marker_expression_plots_noAxis.pdf", sep=""), plot = marker_plot, width=9, height=9 )

# basic QC of data set
VlnPlot(object = seurat.object[[ 2 ]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "genotype")
ggsave( paste(base, "/plots/Supplement/FigS6E_P1_basic_QC.pdf", sep=""), width=7, height=4)

# ========================================================
# extract relevant cells, focusing on excitatory neurons
# ========================================================

Idents ( seurat.object[[ 1 ]] ) <- "seurat_clusters"
seurat.object[[ 1 ]] <- subset( seurat.object[[ 1 ]], idents = c(  "1", "9", "7", "3", "2", "0" ), invert = F )
Idents( seurat.object[[ 1 ]] ) <- "predicted.id"
seurat.object[[ 1 ]] <- subset( seurat.object[[ 1 ]], idents = "Migrating_neurons_P1" )

seurat.object[[ 1 ]]
# An object of class Seurat 
# 20419 features across 5296 samples within 1 assay 
# Active assay: RNA (20419 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

Idents ( seurat.object[[ 2 ]] ) <- "seurat_clusters"
seurat.object[[ 2 ]] <- subset( seurat.object[[ 2 ]], idents = c(  "4", "0" ), invert = F )
Idents( seurat.object[[ 2 ]] ) <- "predicted.id"
seurat.object[[ 2 ]] <- subset( seurat.object[[ 2 ]], idents = "Migrating_neurons_P1" )

seurat.object[[ 2 ]]
# An object of class Seurat 
# 20851 features across 2148 samples within 1 assay 
# Active assay: RNA (20851 features, 2000 variable features)
# 3 layers present: counts, data, scale.data
# 2 dimensional reductions calculated: pca, umap

# redo standard pipeline for filtered dataset
E18P1_Rptor_dataset <- merge( seurat.object[[ 1 ]], seurat.object[[ 2 ]] )

# An object of class Seurat 
# An object of class Seurat 
# 21481 features across 7444 samples within 1 assay 
# Active assay: RNA (21481 features, 2000 variable features)
# 6 layers present: counts.1, counts.2, data.1, scale.data.1, data.2, scale.data.2

# filter for unambiguous hash tags
Idents( E18P1_Rptor_dataset ) <- "hto_demux"
E18P1_Rptor_dataset <- subset( E18P1_Rptor_dataset, idents = c(  "Doublet", "Negative" ), invert = T )

# An object of class Seurat 
# 21481 features across 5241 samples within 1 assay 
# Active assay: RNA (21481 features, 2000 variable features)
# 6 layers present: counts.1, counts.2, data.1, scale.data.1, data.2, scale.data.2

# standard workflow
E18P1_Rptor_dataset <- JoinLayers( E18P1_Rptor_dataset )
E18P1_Rptor_dataset <- FindVariableFeatures( E18P1_Rptor_dataset )
E18P1_Rptor_dataset <- ScaleData( E18P1_Rptor_dataset )
E18P1_Rptor_dataset <- RunPCA( E18P1_Rptor_dataset )

ElbowPlot( E18P1_Rptor_dataset )
ggsave( filename = paste(base, "plots/QC/E18_P1_combined_filtered_elbow.png", sep=""), width = 5, height = 5 )

E18P1_Rptor_dataset <- RunUMAP( E18P1_Rptor_dataset, dims = 1:15 )
DimPlot( object = E18P1_Rptor_dataset, reduction = "umap")
ggsave( filename = paste(base, "plots/QC/E18_P1_combined_filtered_UMAP.png", sep=""), width = 5, height = 5 )

# ==================================================
# classify cells in UL/DL based on marker expression
# ==================================================
reference@meta.data$mod_celltype <- reference@meta.data$New_cellType
# DL cells separate into multiple types, combine here
reference@meta.data$mod_celltype[ which( reference@meta.data$mod_celltype %in% c("CThPN_P1", "SCPN_P1", "Layer_4_P1", "Layer_6b_P1", "NP_P1" ) )] <- "DL_CPN_P1"
Idents(reference) <- "mod_celltype"

DimPlot( reference, reduction = "umap", group.by = "mod_celltype", label=T)
ggsave( paste(base, "plots/QC/Reference_P1_umap_mod_celltype.pdf", sep=""), width=10, height = 7)

layer_markers_ref <- FindMarkers( object = reference, ident.1 = "UL_CPN_P1", ident.2 = "DL_CPN_P1")
layer_markers_ref <- layer_markers_ref[ order(layer_markers_ref$p_val_adj, decreasing = F),]
layer_markers_ref_pos <- layer_markers_ref[ which(layer_markers_ref$avg_log2FC > 0), ]
layer_markers_ref_neg <- layer_markers_ref[ which(layer_markers_ref$avg_log2FC < 0), ]

## classify cells as Ctip2 and Satb2
E18P1_Rptor_dataset <-  AddModuleScore( object = E18P1_Rptor_dataset, features = list( c("Bcl11b", "Foxp2", "Fezf2"), c("Satb2", "Cux1", "Cux2") ) )

#plot ModuleScore cutoffs
ggplot(E18P1_Rptor_dataset@meta.data, aes(x=Cluster1, y=Cluster2)) + geom_point(size=0.5) + 
  geom_hline(yintercept = -0.9) + geom_vline(xintercept = 0) + theme_classic()
ggsave( filename = paste(base, "plots/QC/E18_P1_combined_filtered_cluster_dotplot.png", sep=""), width = 5, height = 5 )

# classification based on cutoffs
E18P1_Rptor_dataset@meta.data$marker_class <- NA
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 > 0 & E18P1_Rptor_dataset@meta.data$Cluster2 < -0.9 ) ] <- "Ctip2"
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 < 0 & E18P1_Rptor_dataset@meta.data$Cluster2 > -0.9 ) ] <- "Satb2"
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 > 0 & E18P1_Rptor_dataset@meta.data$Cluster2 > -0.9 ) ] <- "both"
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 < 0 & E18P1_Rptor_dataset@meta.data$Cluster2 < -0.9 ) ] <- "none"

# ==================================================
# refine classification based on reference markers
# ==================================================

Idents( E18P1_Rptor_dataset ) <- "marker_class"
# define layer markers in E18/P1 data
layer_markers <- FindMarkers( object = E18P1_Rptor_dataset, ident.1 = "Satb2", ident.2 = "Ctip2")
layer_markers <- layer_markers[ order(layer_markers$p_val_adj, decreasing = F),]
layer_markers_pos <- layer_markers[ which(layer_markers$avg_log2FC > 0), ]
layer_markers_neg <- layer_markers[ which(layer_markers$avg_log2FC < 0), ]
# overlap both datasets to find a few more markers for cleaner classification
refined_markers_pos <- intersect( rownames(layer_markers_pos)[1:20], rownames(layer_markers_ref_pos) )
#"Satb2"  "Cux1"   "Ppfia2" "Frmd4a" "Cux2"   "Robo2"  "Gria2"  "Hs6st2" "Clstn2" "Rbfox1" "Ptprz1" "Gpm6a"  "Ptn"    "Ttc28"  "Osbpl6"
refined_markers_neg <- intersect( rownames(layer_markers_neg)[1:20], rownames(layer_markers_ref_neg) )
#"Bcl11b" "Foxp2"  "Fezf2"  "Pcp4"   "Rprm"   "Synpr"  "Zfpm2"  "Rnd2"  

E18P1_Rptor_dataset <-  AddModuleScore( object = E18P1_Rptor_dataset, features = list( refined_markers_pos, refined_markers_neg ) )

#plot ModuleScore cutoffs
ggplot(E18P1_Rptor_dataset@meta.data, aes(x=Cluster1, y=Cluster2)) + geom_point(size=0.5) + 
  geom_hline(yintercept = -0.2) + geom_vline(xintercept = 0) + theme_classic()
ggsave( filename = paste(base, "plots/QC/E18_P1_combined_filtered_cluster_dotplot_refined.png", sep=""), width = 5, height = 5 )

E18P1_Rptor_dataset@meta.data$marker_class <- NA
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 > 0 & E18P1_Rptor_dataset@meta.data$Cluster2 < -0.2 ) ] <- "Satb2"
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 < 0 & E18P1_Rptor_dataset@meta.data$Cluster2 > -0.2 ) ] <- "Ctip2"
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 > 0 & E18P1_Rptor_dataset@meta.data$Cluster2 > -0.2 ) ] <- "both"
E18P1_Rptor_dataset@meta.data$marker_class[ which( E18P1_Rptor_dataset@meta.data$Cluster1 < 0 & E18P1_Rptor_dataset@meta.data$Cluster2 < -0.2 ) ] <- "none"

DimPlot( E18P1_Rptor_dataset, reduction = "umap", group.by = "marker_class", label=F)

ggsave( filename = paste(base, "plots/QC/E18_P1_combined_filtered_cluster_UMAP.png", sep=""), width = 5, height = 5 )

# modify genotype to be a more useful
E18P1_Rptor_dataset@meta.data$genotype <- sapply ( E18P1_Rptor_dataset@meta.data$genotype, function (x) strsplit(x = x, split = "_")[[1]][2] )
# extract those neurons that were assigned unambiguously
Idents( E18P1_Rptor_dataset ) <- "marker_class"
E18P1_Rptor_dataset <- subset( x = E18P1_Rptor_dataset, idents = c("Ctip2", "Satb2") )

DimPlot( E18P1_Rptor_dataset, reduction = "umap", group.by = "marker_class", label=F)
ggsave( filename = paste(base, "plots/QC/E18_P1_neurons_pure_UMAP.png", sep=""), width = 5, height = 4 )

# An object of class Seurat 
# 21481 features across 2819 samples within 1 assay 
# Active assay: RNA (21481 features, 2000 variable features)
# 3 layers present: scale.data, data, counts
# 2 dimensional reductions calculated: pca, umap

# define layer markers in E18/P1 data
layer_markers <- FindMarkers( object = E18P1_Rptor_dataset, ident.1 = "Satb2", ident.2 = "Ctip2")
tmp_gl <- layer_markers[which(layer_markers$p_val_adj < 1e-20),]
tmp_gl <- tmp_gl[ which( rownames (tmp_gl) %in% rownames( E18P1_Rptor_dataset@assays$RNA$scale.data ) ),]
marker_gl <- rownames(tmp_gl[order(tmp_gl$avg_log2FC, decreasing = T),])

marker_heatmap <- DoHeatmap( object = E18P1_Rptor_dataset, features = marker_gl)
ggsave( plot = marker_heatmap, filename = paste(base, "plots/QC/E18_P1_neuron_classes_marker_heatmap.pdf", sep=""), width = 7, height = 13 )

idx2keep <- which( marker_gl %in% c("Bcl11b", "Foxp2", "Fezf2", "Satb2", "Cux1", "Cux2") )
new_label <- rep("", length(marker_gl))
new_label[ idx2keep ] <- marker_gl[idx2keep]

marker_heatmap <- marker_heatmap + scale_y_discrete(label=rev(new_label) )
ggsave( plot = marker_heatmap, filename = paste(base, "plots/Supplement/FigS6I_E18_P1_neuron_classes_marker_heatmap.pdf", sep=""), width = 7, height = 13 )

# plot relative distribution of cells per marker class and genotype
cdata <- ddply(E18P1_Rptor_dataset@meta.data, c("age", "genotype", "marker_class"), summarise,
               N    = length(marker_class)
)

cdata$age_gen <- paste( cdata$age, cdata$genotype, sep="_" )

E18P1_Rptor_dataset@meta.data$age_gen <- paste( E18P1_Rptor_dataset@meta.data$age, E18P1_Rptor_dataset@meta.data$genotype, sep="_")
tot_neurons <- table( E18P1_Rptor_dataset$age_gen)

cdata$tot <-tot_neurons[ cdata$age_gen ]
cdata$rel <- cdata$N / cdata$tot

ggplot( cdata, aes(x=marker_class, y=rel)) + 
  geom_bar(stat="identity", position = "dodge") + facet_grid(~age_gen) + theme_classic()
ggsave( filename = paste(base, "plots/QC/E18_P1_neuron_classes_rel_count_age_gen.pdf", sep=""), width = 7, height = 3 )

ggplot( cdata, aes(x=marker_class, y=N)) + 
  geom_bar(stat="identity", position = "dodge") + facet_grid(~age_gen) + theme_classic()
ggsave( filename = paste(base, "plots/QC/E18_P1_neuron_classes_abs_count_age_gen.pdf", sep=""), width = 7, height = 3 )

# plot layer abundance independent of age
cdata <- ddply(E18P1_Rptor_dataset@meta.data, c("genotype", "marker_class"), summarise,
               N    = length(marker_class)
)

tot_neurons <- table( E18P1_Rptor_dataset$genotype)

cdata$tot <-tot_neurons[ cdata$genotype ]
cdata$rel <- cdata$N / cdata$tot
cdata$marker_class <- factor( cdata$marker_class, levels = c("Satb2", "Ctip2"))

conf_int <- t(sapply (1:nrow(cdata), function (x) {
  tmp <- clopper.pearson.ci(cdata$N[x], cdata$tot[x], alpha = 0.05, CI = "two.sided")
}))
conf_int <- as.data.frame( conf_int )

cdata <- cbind( cdata, conf_int )
cdata$Lower.limit <- as.numeric( cdata$Lower.limit )
cdata$Upper.limit <- as.numeric( cdata$Upper.limit )

ggplot( cdata, aes(x=marker_class, y=rel)) + 
  geom_bar(stat="identity", position = "dodge") + geom_errorbar( aes(ymin = Lower.limit, ymax = Upper.limit)) +
  facet_grid(~genotype) + theme_classic()
ggsave( filename = paste(base, "plots/Supplement/FigS6J_E18_P1_neuron_classes_gen_rel.pdf", sep=""), width = 5, height = 5 )

E18P1_Rptor_dataset@meta.data$age_markerClass <- paste( E18P1_Rptor_dataset$age, E18P1_Rptor_dataset$marker_class, sep="_" )

# pseudobulk analysis to plot marker gene expression
pb_mat <- AggregateExpression(
  E18P1_Rptor_dataset,
  group.by = "age_markerClass", return.seurat = T)

library(pheatmap)
pheatmap( pb_mat@assays$RNA$scale.data[c("Bcl11b", "Foxp2", "Fezf2", "Pcp4", "Satb2", "Cux1", "Cux2", "Robo2"),],
          filename =  paste(base, "plots/Main/Fig6C_neuron_layer_marker_expression.pdf", sep="") )
# sometimes necessary after pheatmap call
dev.off()

# ==============
# DEG analysis
# ==============

# DEGs between genotypes in each layer, separately for E18 and P1
group_deg_list <- list()
for ( age in unique( E18P1_Rptor_dataset@meta.data$age ) ) {
  
  group_deg_list[[ age ]] <- list()
  
  for ( layer in unique( E18P1_Rptor_dataset@meta.data$marker_class ) ) {
    
    Idents( E18P1_Rptor_dataset ) <- "marker_class"
    tmp_obj <- subset( x = E18P1_Rptor_dataset, idents = layer )
    Idents( tmp_obj ) <- "age"
    tmp_obj <- subset( x = tmp_obj, idents = age )
    Idents( tmp_obj ) <- "genotype"
    group_deg_list[[ age ]][[ layer]] <- FindMarkers( tmp_obj, ident.1 = "Het", ident.2 = "WT", test.use = "MAST", logfc.threshold = 0.01)
    
  }
  
}

saveRDS( object = group_deg_list, file = paste(base, "RDS_files/E18_P1_group_deg_list.RDS", sep=""))

# DEGS between upper layer and lower layer separately for genotype and age
layer_deg_list <- list()
for ( age in unique( E18P1_Rptor_dataset@meta.data$age ) ) {
  
  layer_deg_list[[ age ]] <- list()
  
  for ( genotype in unique( E18P1_Rptor_dataset@meta.data$genotype ) ) {
    
    Idents( E18P1_Rptor_dataset ) <- "genotype"
    tmp_obj <- subset( x = E18P1_Rptor_dataset, idents = genotype )
    Idents( tmp_obj ) <- "age"
    tmp_obj <- subset( x = tmp_obj, idents = age )
    Idents( tmp_obj ) <- "marker_class"
    layer_deg_list[[ age ]][[ genotype]] <- FindMarkers( tmp_obj, ident.1 = "Satb2", ident.2 = "Ctip2", test.use = "MAST", logfc.threshold = 0.01)
    
  }
  
}

saveRDS( object = layer_deg_list, file = paste(base, "RDS_files/E18_P1_layer_deg_list.RDS", sep=""))

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
#   [1] pheatmap_1.0.12        GeneOverlap_1.38.0     ggplot2_3.5.0          plyr_1.8.9             Seurat_5.0.1          
# [6] SeuratObject_5.0.1     sp_2.1-3               org.Mm.eg.db_3.18.0    AnnotationDbi_1.64.1   IRanges_2.36.0        
# [11] S4Vectors_0.40.2       Biobase_2.62.0         BiocGenerics_0.48.1    clusterProfiler_4.10.1
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.3                    matrixStats_1.2.0           spatstat.sparse_3.0-3       bitops_1.0-7               
# [5] enrichplot_1.22.0           HDO.db_0.99.1               httr_1.4.7                  RColorBrewer_1.1-3         
# [9] tools_4.3.2                 sctransform_0.4.1           utf8_1.2.4                  R6_2.5.1                   
# [13] lazyeval_0.2.2              uwot_0.1.16                 withr_3.0.0                 prettyunits_1.2.0          
# [17] gridExtra_2.3               progressr_0.14.0            cli_3.6.2                   textshaping_0.3.7          
# [21] spatstat.explore_3.2-6      fastDummies_1.7.3           scatterpie_0.2.2            labeling_0.4.3             
# [25] spatstat.data_3.0-4         ggridges_0.5.6              pbapply_1.7-2               systemfonts_1.0.6          
# [29] yulab.utils_0.1.4           gson_0.1.0                  DOSE_3.28.2                 parallelly_1.37.0          
# [33] limma_3.58.1                rstudioapi_0.16.0           RSQLite_2.3.6               generics_0.1.3             
# [37] gridGraphics_0.5-1          gtools_3.9.5                ica_1.0-3                   spatstat.random_3.2-2      
# [41] dplyr_1.1.4                 GO.db_3.18.0                Matrix_1.6-5                ggbeeswarm_0.7.2           
# [45] fansi_1.0.6                 abind_1.4-5                 lifecycle_1.0.4             SummarizedExperiment_1.32.0
# [49] gplots_3.1.3.1              qvalue_2.34.0               SparseArray_1.2.4           Rtsne_0.17                 
# [53] grid_4.3.2                  blob_1.2.4                  promises_1.2.1              crayon_1.5.2               
# [57] miniUI_0.1.1.1              lattice_0.21-9              cowplot_1.1.3               KEGGREST_1.42.0            
# [61] pillar_1.9.0                fgsea_1.28.0                GenomicRanges_1.54.1        future.apply_1.11.1        
# [65] codetools_0.2-19            fastmatch_1.1-4             leiden_0.4.3.1              glue_1.7.0                 
# [69] packrat_0.9.2               ggfun_0.1.4                 data.table_1.15.0           vctrs_0.6.5                
# [73] png_0.1-8                   treeio_1.26.0               spam_2.10-0                 gtable_0.3.4               
# [77] cachem_1.0.8                S4Arrays_1.2.1              mime_0.12                   tidygraph_1.3.1            
# [81] survival_3.5-7              SingleCellExperiment_1.24.0 statmod_1.5.0               ellipsis_0.3.2             
# [85] fitdistrplus_1.1-11         ROCR_1.0-11                 nlme_3.1-163                ggtree_3.10.1              
# [89] bit64_4.0.5                 progress_1.2.3              RcppAnnoy_0.0.22            GenomeInfoDb_1.38.8        
# [93] irlba_2.3.5.1               vipor_0.4.7                 KernSmooth_2.23-22          colorspace_2.1-0           
# [97] DBI_1.2.2                   ggrastr_1.0.2               tidyselect_1.2.1            bit_4.0.5                  
# [101] compiler_4.3.2              hdf5r_1.3.10                DelayedArray_0.28.0         plotly_4.10.4              
# [105] shadowtext_0.1.3            scales_1.3.0                caTools_1.18.2              lmtest_0.9-40              
# [109] stringr_1.5.1               digest_0.6.34               goftest_1.2-3               spatstat.utils_3.1-0       
# [113] presto_1.0.0                XVector_0.42.0              htmltools_0.5.7             pkgconfig_2.0.3            
# [117] MatrixGenerics_1.14.0       fastmap_1.1.1               rlang_1.1.3                 htmlwidgets_1.6.4          
# [121] shiny_1.8.0                 farver_2.1.1                zoo_1.8-12                  jsonlite_1.8.8             
# [125] BiocParallel_1.36.0         GOSemSim_2.28.1             RCurl_1.98-1.14             magrittr_2.0.3             
# [129] GenomeInfoDbData_1.2.11     ggplotify_0.1.2             dotCall64_1.1-1             patchwork_1.2.0            
# [133] munsell_0.5.0               Rcpp_1.0.12                 ape_5.8                     viridis_0.6.5              
# [137] reticulate_1.35.0           stringi_1.8.3               ggraph_2.2.1                zlibbioc_1.48.2            
# [141] MASS_7.3-60                 MAST_1.28.0                 parallel_4.3.2              listenv_0.9.1              
# [145] ggrepel_0.9.5               deldir_2.0-2                Biostrings_2.70.3           graphlayouts_1.1.1         
# [149] splines_4.3.2               tensor_1.5                  hms_1.1.3                   igraph_2.1.4               
# [153] spatstat.geom_3.2-8         RcppHNSW_0.6.0              reshape2_1.4.4              tweenr_2.0.3               
# [157] httpuv_1.6.14               RANN_2.6.1                  tidyr_1.3.1                 purrr_1.0.2                
# [161] polyclip_1.10-6             future_1.33.1               scattermore_1.2             ggforce_0.4.2              
# [165] xtable_1.8-4                RSpectra_0.16-1             tidytree_0.4.6              later_1.3.2                
# [169] viridisLite_0.4.2           ragg_1.3.0                  tibble_3.2.1                aplot_0.2.2                
# [173] memoise_2.0.1               beeswarm_0.4.0              cluster_2.1.4               globals_0.16.2  
