# Import packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)

# Configuration
set.seed(1234)
setwd("/Users/nguyenpu/Documents/GitHub/serezani_lab/10x_visium/")
datadir <- "results/010_Normalization/"
outputdir <- "results/020_Clustering/"

if (!dir.exists(outputdir)) {
    dir.create(outputdir)
}

# Upload normalized counts
ab1 <- readRDS(paste0(datadir, "ab1_normalized.rds"))
ab2 <- readRDS(paste0(datadir, "ab2_normalized.rds"))
ab3 <- readRDS(paste0(datadir, "ab3_normalized.rds"))
ab4 <- readRDS(paste0(datadir, "ab4_normalized.rds"))
ab5 <- readRDS(paste0(datadir, "ab5_normalized.rds"))
ab6 <- readRDS(paste0(datadir, "ab6_normalized.rds"))
ab7 <- readRDS(paste0(datadir, "ab7_normalized.rds"))
ab8 <- readRDS(paste0(datadir, "ab8_normalized.rds"))

# Store as a list
abs <- list(ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8)

###### Dimensionality Reduction #####
dimred <- function(ab){
    ab <- RunPCA(ab, assay = "SCT", verbose = FALSE)
    ab <- FindNeighbors(ab, reduction = "pca", dims = 1:30)
    ab <- FindClusters(ab, resolution = 0.6, verbose = FALSE) #default resolution 0.8, change as needed
    ab <- RunUMAP(ab, reduction = "pca", dims = 1:30)
    return(ab) 
}

# Test-Run with AB1
ab1 <- dimred(ab1)
p1 <- DimPlot(ab1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ab1, label = TRUE, label.size = 3, pt.size.factor = 1.6)
p1 + p2

# Loop through all samples
for (i in seq_along(abs)) {
    print(paste0("Starting dimensional analysis for AB", i, "."))
    abs[[i]] <- dimred(abs[[i]])  
    saveRDS(abs[[i]], file = paste0(outputdir, "Resolution_0.6/", "AB", i, "_DR.rds"))
    print(paste0("Dimensional analysis for AB", i, " completed."))

    # Plot and Save
    print("Now plotting and saving")
    p1 <- DimPlot(abs[[i]], reduction = "umap", label = TRUE)
    p2 <- SpatialDimPlot(abs[[i]], label = TRUE, label.size = 3, pt.size.factor = 1.6)
    ggsave(paste0(outputdir, "Resolution_0.6/", "UMAP+spatial_plot-AB", i, ".png"),
           plot = p1 + p2,
           units = "in", width = 8, height = 4, scale = 1, dpi = 300)
}

# Interactive plots
SpatialDimPlot(ab1, interactive = TRUE)
ab1 <- SetAllIndent(ab1, id)
