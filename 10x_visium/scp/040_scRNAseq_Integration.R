# Load packages
library(SeuratDisk)
library(Seurat)
library(zellkonverter)
library(scRNAseq)
library(scran)
library(scater)
library(patchwork)
library(dplyr)

# Configuration
set.seed(1234)
setwd("/Volumes/SEREZANI/001_Phoebe/10x_visium/")

datadir <- "counts/skin_scrna_seq_references/"
outputdir <- "results/040_scRNAseq_Integration/"

if (dir.exists(outputdir) = FALSE){
    dir.create(outputdir)
}


##### Almat et al #####

# Load nondiabetic integrated skin datasets for cell-typing
# Source: https://zenodo.org/records/10198892
# Convert the scanpy .h5ad file to a single cell experiment
sce <- readH5AD(paste0(datadir, "integratedskindata.h5ad"))
saveRDS(sce, file = paste0(outputdir, "integrated_skin_reference-sce.rds"))

# Convert to Seurat, specifying the correct assays
integrated_skin <- as.Seurat(sce, counts="counts", data="X")
# Rename the assay "originalexp" to "RNA"
integrated_skin <- RenameAssays(integrated_skin, originalexp = "RNA")


# saveRDS(integrated_skin, file = paste0(outputdir, "integrated_skin_reference.rds"))
#Before merging common cell types together
p <- DimPlot(integrated_skin, reduction = "X_umap", group.by = "leiden", label = TRUE)
p

#Merge common cell types/leidin clusters together
table(integrated_skin$leiden)

# Ensure leiden is a character vector before renaming
integrated_skin$leiden_merged <- as.character(integrated_skin$leiden)

# Define the renaming mapping
rename_mapping <- list("EPI-1" = "Epidermal",
                       "EPI-2" = "Epidermal",
                       "EPI-3" = "Epidermal",
                       "EPI-4" = "Epidermal",
                       "EPI-5" = "Epidermal",
                       "FIB-1" = "Fibroblasts",
                       "FIB-2" = "Fibroblasts",
                       "FIB-3" = "Fibroblasts",
                       "FIB-4" = "Fibroblasts",
                       "IMM-1" = "Immune",
                       "IMM-2" = "Immune",
                       "IMM-3" = "Immune",
                       "PERI" = "Pericytes",
                       "ENDO" = "Endothelial",
                       "SMC" = "Smooth Muscle Cells",
                       "LYM-E" = "Lymphatic Endothelial",
                       "SCH" = "Schwann Cells",
                       "MEL" = "Melanocytes",
                       "LC" = "Langerhans Cells",
                       "SKEL" = "Skeletal Muscle")

# Apply renaming safely
integrated_skin$leiden_merged <- sapply(integrated_skin$leiden_merged, function(x) {
  if (x %in% names(rename_mapping)) {
    rename_mapping[[x]]
  } else {
    x  # Keep original name if not in mapping
  }
})

# Convert back to factor (if needed)
integrated_skin$leiden_merged <- factor(integrated_skin$leiden_merged)

# Verify the renaming
table(integrated_skin$leiden)

#After merging common cell types together
p <- DimPlot(integrated_skin, reduction = "X_umap", group.by = "leiden", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(paste0(outputdir, "integrated_skin-UMAP.png"), plot = p, height = 10, width = 10, units = "in", dpi = 300)

saveRDS(integrated_skin, file = paste0(outputdir, "integrated_skin_reference-Seurat.rds"))

integrated_skin <- readRDS(paste0(outputdir, "integrated_skin_reference-Seurat.rds"))
# SCT Transform the integrated_skin reference
# Not enough RAM, try another dataset
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4GB
integrated_skin <- SCTransform(integrated_skin, n_cells = 3000, verbose = FALSE)


##### Guerrero-Juarez et al #####
# https://pubmed.ncbi.nlm.nih.gov/30737373/
guerrero_juarez <- CreateSeuratObject(Read10X(paste0(datadir, "Guerrero-Juarez_2019/GSE113854")))
saveRDS(guerrero_juarez, file = paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))
guerrero_juarez <- readRDS(paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))

# QC
guerrero_juarez[["percent.mt"]] <- PercentageFeatureSet(guerrero_juarez, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(guerrero_juarez, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0(outputdir, "Guerrero_Juarez-2019/guerrero_juarez-QC.png"),
        plot = last_plot(),
        height = 10, width = 10,
        units = "in",
        dpi = 300)

#  Quality control metrics included keeping cells displaying < 8000 UMI/cell and < 2500 genes/cell, 
#  and no more than 8% mitochondrial gene expression

guerrero_juarez <- subset(guerrero_juarez,
                          subset = nCount_RNA < 8000 &
                          nFeature_RNA < 2500 &
                          percent.mt < 8)
dim(guerrero_juarez)

# SCT Transform the guerrero_juarez data
options(future.globals.maxSize = 4 * 1024^3)
guerrero_juarez <- SCTransform(guerrero_juarez, ncells = 3000, verbose = FALSE)
# guerrero_juarez <- NormalizeData(guerrero_juarez, verbose = FALSE)

# Save the processed data
saveRDS(guerrero_juarez, file = paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))

# Repeat dimensional reduction for cell clustering but with SCTransform instead of normal Log2Norm
guerrero_juarez <- FindVariableFeatures(guerrero_juarez, selection.method = "vst", nfeatures = 2000)
guerrero_juarez <- ScaleData(guerrero_juarez, features = rownames(guerrero_juarez))
guerrero_juarez <- RunPCA(guerrero_juarez, verbose = FALSE)

# Jackstraw doesn't work on SCT transfromed Data
# guerrero_juarez <- JackStraw(guerrero_juarez, num.replicate = 100, dims = 50)
# guerrero_juarez <- ScoreJackStraw(guerrero_juarez, dims = 1:50)

# saveRDS(guerrero_juarez, file = paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))
# Load guerrero_juarez Seurat object from 03/13
# guerrero_juarez <- readRDS(paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))


ElbowPlot(guerrero_juarez, ndims = 40) 
guerrero_juarez <- FindNeighbors(guerrero_juarez, dims = 1:30, verbose = FALSE)
guerrero_juarez <- FindClusters(guerrero_juarez, resolution = 0.35, verbose = FALSE)

guerrero_juarez <- RunUMAP(guerrero_juarez, dims = 1:30)
guerrero_juarez <- RunTSNE(guerrero_juarez, dims = 1:30)

DimPlot(guerrero_juarez, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(guerrero_juarez, reduction = "tsne", label = TRUE, pt.size = 0.5)
ggsave(filename = paste0(outputdir, "Guerrero_Juarez-2019/TSNE_Plot_unlabeled.png"),
       width = 10,
       height = 10,  # Changed "length" to "height"
       units = "in")

saveRDS(guerrero_juarez, file = paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))

# Find cluster markers and compare to reference

cluster_markers <- FindAllMarkers(guerrero_juarez, 
                                  only.pos = TRUE, 
                                  min.pct = 0.1, 
                                  logfc.threshold = 0.25,
                                  test.use = "LR")  # Likelihood Ratio Test

cluster_markers <- cluster_markers[cluster_markers$p_val < 0.01, ]
write.csv(cluster_markers, file = paste0(outputdir, "Guerrero_Juarez-2019/TSNE_ClusterMarkers-SCT.csv"))

guerrero_juarez <- readRDS(paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))
cluster_markers <- read.csv(paste0(outputdir, "Guerrero_Juarez-2019/TSNE_ClusterMarkers-SCT.csv"), header = TRUE)

# Rename clusters:
# Define the Leiden-to-cell type mapping
# Define the mapping of clusters based on Guerrero-Juarez et al. (2019)
cluster_mapping <- c(
    "0" = "Fibroblasts, type 1 (FIB-1)", "1" = "Fibroblasts, type 4 (FIB-4)",
    "2" = "Myeloid Cells (MYL)", "3" = "Fibroblasts, type 3 (FIB-3)",
    "4" = "Endothelial Cells (ENDO)", "5" = "Dendritic Cells (DEN)",
    "6" = "Fibroblasts, type 1 (FIB-1)", "7" = "Dendritic Cells (DEN)",
    "8" = "Dendritic Cells (DEN)", "9" = "B Lymphocytes (BCELL)",
    "10" = "Fibroblasts, type 3 (FIB-3)", "11" = "B Lymphocytes (BCELL)",
    "12" = "Schwann Cells (SCH)", "13" = "Dendritic Cells (DEN)",
    "14" = "Lymphatic Endothelial Cells (LYME)"
)

# Rename clusters in the Seurat object
guerrero_juarez$renamed_cluster <- sapply(guerrero_juarez$seurat_clusters, function(x) {
  if (x %in% names(cluster_mapping)) {
    cluster_mapping[[as.character(x)]]
  } else {
    x  # Keep original name if not in mapping
  }
})

# Verify the renaming
table(guerrero_juarez$renamed_cluster)

#After merging common cell types together
p <- DimPlot(guerrero_juarez, reduction = "umap", group.by = "renamed_cluster", label = TRUE, pt.size = 0.5) + NoLegend()
p
saveRDS(guerrero_juarez, file = paste0(datadir, "Guerrero-Juarez_2019/GSE113854/GSE113854_SeuratObj.rds"))



##### Integration with scRNAseq data #####

refdir <- ("counts/skin_scrna_seq_references/Guerrero-Juarez_2019/")
datadir <- ("results/020_Clustering/")

#Load in Guerrero-Juarez reference
guerrero_juarez <- readRDS(paste0(refdir, "GSE113854/GSE113854_SeuratObj.rds"))

# Load in the clustered 10X Visium Data
ab1 <- readRDS(paste0(datadir, "AB1_DR.rds"))
ab2 <- readRDS(paste0(datadir, "AB2_DR.rds"))
ab3 <- readRDS(paste0(datadir, "AB3_DR.rds"))
ab4 <- readRDS(paste0(datadir, "AB4_DR.rds"))
ab5 <- readRDS(paste0(datadir, "AB5_DR.rds"))
ab6 <- readRDS(paste0(datadir, "AB6_DR.rds"))
ab7 <- readRDS(paste0(datadir, "AB7_DR.rds"))
ab8 <- readRDS(paste0(datadir, "AB8_DR.rds"))

common_features <- intersect(rownames(guerrero_juarez), rownames(ab1))
sct_features_ref <- rownames(guerrero_juarez[["SCT"]]@scale.data)
sct_features_query <- rownames(ab1[["SCT"]]@scale.data)

# Find intersection of common features and scaled SCT features
final_features <- intersect(common_features, intersect(sct_features_ref, sct_features_query))
length(final_features)  # See how many are left

seurat_list <- PrepSCTIntegration(
   object.list = list(guerrero_juarez, ab1), 
   anchor.features = final_features
)

# Ensure ab1 remains a Seurat object
ab1 <- seurat_list[[2]]  # Extract the modified query object

anchors <- FindTransferAnchors(
   reference = guerrero_juarez,
   query = ab1,
   normalization.method = "SCT",
   reference.assay = "SCT",
   query.assay = "SCT",
   features = NULL,
)

predictions.assay <- TransferData(
   anchorset = anchors, 
   refdata = guerrero_juarez$renamed_cluster,
   prediction.assay = TRUE,
   weight.reduction = guerrero_juarez[["pca"]],
   dims = 1:min(30, ncol(Embeddings(guerrero_juarez[["pca"]])))
)

