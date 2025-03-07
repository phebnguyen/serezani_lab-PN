#The purpose of this script is to normalize the spatial data post-alignment by CellRanger

# Import packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Configuration
set.seed(1234)
setwd("/Users/nguyenpu/Documents/GitHub/serezani_lab/10x_visium/")
datadir <- "counts/"
outputdir <- "results/010_Normalization/"

if (dir.exists(outputdir) != TRUE) {
    dir.create(outputdir)
}

# Load the spatial data
# Key: AB1-4 = Non-Diabetic
ab1 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-1/"))
ab2 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-2/"))
ab3 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-3/"))
ab4 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-4/"))
# Key: AB5-8 = Diabetic
ab5 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-5/"))
ab6 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-6/"))
ab7 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-7/"))
ab8 <- Load10X_Spatial(data.dir = paste0(datadir, "count_10516-AB-8/"))

# Normalize the spatial data
ab1 <- SCTransform(ab1, assay = "Spatial", verbose = FALSE)
ab2 <- SCTransform(ab2, assay = "Spatial", verbose = FALSE)
ab3 <- SCTransform(ab3, assay = "Spatial", verbose = FALSE)
ab4 <- SCTransform(ab4, assay = "Spatial", verbose = FALSE)
ab5 <- SCTransform(ab5, assay = "Spatial", verbose = FALSE)
ab6 <- SCTransform(ab6, assay = "Spatial", verbose = FALSE)
ab7 <- SCTransform(ab7, assay = "Spatial", verbose = FALSE)
ab8 <- SCTransform(ab8, assay = "Spatial", verbose = FALSE)

# Save the normalized data to an RDS file
saveRDS(ab1, file = paste0(outputdir, "ab1_normalized.rds"))
saveRDS(ab2, file = paste0(outputdir, "ab2_normalized.rds"))
saveRDS(ab3, file = paste0(outputdir, "ab3_normalized.rds"))
saveRDS(ab4, file = paste0(outputdir, "ab4_normalized.rds"))
saveRDS(ab5, file = paste0(outputdir, "ab5_normalized.rds"))
saveRDS(ab6, file = paste0(outputdir, "ab6_normalized.rds"))
saveRDS(ab7, file = paste0(outputdir, "ab7_normalized.rds"))
saveRDS(ab8, file = paste0(outputdir, "ab8_normalized.rds"))

# Sanity Check
# Make list of spatial objects and save RDS
abs <- c(ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8)
n = 1
for (ab in abs){
    p <- SpatialFeaturePlot(ab, features = c("Cd68", "Lamp1"), pt.size.factor = 2.6, alpha = c(0.1, 1))
    ggsave(paste0(outputdir,"/AB",n,".png"),
    width = 10, height = 10, scale = 1, units = "in", bg = "white")
    n <- n + 1
}