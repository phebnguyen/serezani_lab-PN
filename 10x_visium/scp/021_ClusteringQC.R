# Purpose: select for the ideal number of clusters to maintain distinct tissue regions that are biologically relevant.

# Import packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(ggtree)

# Configuration
set.seed(1234)
setwd("/Users/nguyenpu/Documents/GitHub/serezani_lab/10x_visium/")
datadir <- "results/020_Clustering/"
outputdir <- "results/020_Clustering/"

# Load dimensionally reduced data
ab1 <- readRDS(paste0(datadir, "AB1_DR.rds"))
ab2 <- readRDS(paste0(datadir, "AB2_DR.rds"))
ab3 <- readRDS(paste0(datadir, "AB3_DR.rds"))
ab4 <- readRDS(paste0(datadir, "AB4_DR.rds"))
ab5 <- readRDS(paste0(datadir, "AB5_DR.rds"))
ab6 <- readRDS(paste0(datadir, "AB6_DR.rds"))
ab7 <- readRDS(paste0(datadir, "AB7_DR.rds"))
ab8 <- readRDS(paste0(datadir, "AB8_DR.rds"))

# Store as a list
abs <- list(ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8)

##### Cluster Trees #####
# Function to generate cluster trees
cluster_tree <- function(abs, resolution){
    resolution <- as.numeric(resolution)

    for (i in seq_along(abs)) {
        print(paste0("Starting cluster tree for AB", i, " with resolution ", resolution, "."))
        ab <- abs[[i]]

        # Set cluster resolution
        res_col <- paste0("SCT_snn_res.", resolution)

        # Set cluster identity to chosen resolution
        Idents(ab) <- res_col

        # Build the cluster tree
        ab <- BuildClusterTree(ab)
        abs[[i]] <- ab  # Save modified object

        # Ensure output directory exists
        res_dir <- paste0(outputdir, "Resolution_", resolution, "/")
        if (!dir.exists(res_dir)) {
            dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
        }

        # Generate the cluster tree plot
        tree <- Tool(ab, slot = "BuildClusterTree")

        png(paste0(res_dir, "ClusterTree_AB", i, ".png"), width = 1000, height = 1250, res = 300)
        plot(tree, main = paste0("Cluster Tree for AB", i, "\n Resolution: ", resolution))
        dev.off()

        print(paste0("Cluster tree for AB", i, " completed."))
    }
}

# Run cluster_tree for multiple resolutions
for (res in c(0.6, 0.8, 1)) {
    cluster_tree(abs, res)
}

##### Spatial Marker Identification #####
# Identification of spatial markers for clusters via differential gene expression against other clusters

# Function to find all cluster_markers
find_cluster_markers <- function(abs, resolution, outputdir) {
    res_col <- paste0("SCT_snn_res.", resolution)

    for (i in seq_along(abs)) {
        ab <- abs[[i]]
        #Set Idents to desired resolution
        Idents(ab) <- res_col
        #Find the number of cluster markers, starting from 0
        cluster_vec <- 0:(length(unique(Idents(ab))) - 1)
        #Iterate through the clusters and find markers for each cluster
        for (cluster in cluster_vec) {
            print(paste0("Finding markers for cluster ", cluster, " in AB", i, "."))
            cluster_markers <- FindMarkers(ab, ident.1 = cluster, ident.2 = NULL, test.use = "wilcox")
            
            # Sort by log2FoldChange (avg_log2FC) in descending order
            cluster_markers <- cluster_markers[order(-cluster_markers$avg_log2FC), ]

            # Save the markers
            write.csv(cluster_markers, file = paste0(outputdir, "Resolution_", resolution,
                                                     "/AB", i, "/Cluster", cluster, "_Markers-", "AB", i, ".csv"))
            print(paste0("Cluster ", cluster, " markers for AB", i, " found." ))
        }
    }
}

for (res in c(0.6, 0.8, 1)) {
    find_cluster_markers(abs, res, outputdir)
}

# # Sort folders in Resolution folders into respective AB folders
# library(stringr)
# for (res in c(0.6, 0.8, 1)) {
#     res_dir <- paste0(outputdir, "Resolution_", res, "/")
#     files <- list.files(res_dir)
#     for (file in files) {
#         ab_num <- str_extract(file, "AB[0-9]+")
#         if (!dir.exists(paste0(res_dir, ab_num))) {
#             dir.create(paste0(res_dir, ab_num), showWarnings = FALSE)
#         }
#         file.rename(from = paste0(res_dir, file), to = paste0(res_dir, ab_num, "/", file))
#     }
# }
