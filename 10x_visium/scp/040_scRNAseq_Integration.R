# Load packages
library(SeuratDisk)
library(Seurat)
library(zellkonverter)
library(scRNAseq)
library(scran)
library(scater)

# Configuration
set.seed(1234)
setwd("/Volumes/SEREZANI/001_Phoebe/10x_visium/")

datadir <- "counts/skin_scrna_seq_references/"
outputdir <- "results/040_scRNAseq_Integration/"

if (dir.exists(outputdir) = FALSE){
    dir.create(outputdir)
}

# Load nondiabetic integrated skin datasets for cell-typing
# Source: https://zenodo.org/records/10198892

# Convert the scanpy .h5ad file to a single cell experiment
sce <- readH5AD(paste0(datadir, "integratedskindata.h5ad"))
sce <- logNormCounts(sce)

# Convert to Seurat, specifying the correct assays
integrated_skin <- as.Seurat(sce)
saveRDS(integrated_skin, file = paste0(outputdir, "integrated_skin_reference.rds"))

#Check clusters for integrated_skin
p <- DimPlot(integrated_skin, reduction = "X_umap", group.by = "leiden", label = TRUE)
ggsave(paste0(outputdir, "integrated_skin-UMAP.png"), plot = p, height = 10, width = 10, units = "in", dpi = 300)

##### Integration with scRNAseq data #####

#Might want to rename leiden to cell type from the integrated_skin reference

anchors <- FindTransferAnchors(reference = integrated_skin, query = ab1, normalization_method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = integrated_skin$leiden, prediction.assay = "TRUE",
                    weight_reduction = ab1[["pca"]], dims = 1:30)