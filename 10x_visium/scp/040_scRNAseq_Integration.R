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
integrated_skin <- as.Seurat(sce, counts="counts", data="logcounts")
# saveRDS(integrated_skin, file = paste0(outputdir, "integrated_skin_reference.rds"))
#Before merging common cell types together
p <- DimPlot(integrated_skin, reduction = "X_umap", group.by = "leiden", label = TRUE)
p

#Merge common cell types/leidin clusters together
table(integrated_skin$leiden)
rename_mapping <- c("EPI-1" = "Epidermal",
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

integrated_skin <- RenameIdents(integrated_skin, rename_mapping)

# Verify the renaming
table(integrated_skin$leiden)

#After merging common cell types together
p <- DimPlot(integrated_skin, reduction = "X_umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(paste0(outputdir, "integrated_skin-UMAP.png"), plot = p, height = 10, width = 10, units = "in", dpi = 300)

saveRDS(integrated_skin, file = paste0(outputdir, "integrated_skin_reference-merged.rds"))

##### Integration with scRNAseq data #####

datadir <- ("results/020_Clustering/")

# Load in the clustered 10X Visium Data
ab1 <- readRDS(paste0(datadir, "AB1_DR.rds"))
ab2 <- readRDS(paste0(datadir, "AB2_DR.rds"))
ab3 <- readRDS(paste0(datadir, "AB3_DR.rds"))
ab4 <- readRDS(paste0(datadir, "AB4_DR.rds"))
ab5 <- readRDS(paste0(datadir, "AB5_DR.rds"))
ab6 <- readRDS(paste0(datadir, "AB6_DR.rds"))
ab7 <- readRDS(paste0(datadir, "AB7_DR.rds"))
ab8 <- readRDS(paste0(datadir, "AB8_DR.rds"))

#Test with AB1
anchors <- FindTransferAnchors(reference = integrated_skin, query = ab1, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = Idents(integrated_skin),
                                  prediction.assay = "TRUE", weight_reduction = ab1[["pca"]], dims = 1:30)