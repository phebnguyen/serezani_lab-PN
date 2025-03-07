# Purpose: subset out UMAP clusters into anatomical regions (e.g., abcess vs. surrounding abcess)

# Import packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)

# Configuration
set.seed(1234)
setwd("/Users/nguyenpu/Documents/GitHub/serezani_lab/10x_visium/")
datadir <- "results/020_Clustering/"
outputdir <- "results/030_AnnatomicalAnnotation/"

if (dir.exists(outputdir)!=TRUE){
    dir.create(outputdir)
}

# Load .rds files from 020_Clustering
# Load dimensionally reduced data
ab1 <- readRDS(paste0(datadir, "AB1_DR.rds"))
ab2 <- readRDS(paste0(datadir, "AB2_DR.rds"))
ab3 <- readRDS(paste0(datadir, "AB3_DR.rds"))
ab4 <- readRDS(paste0(datadir, "AB4_DR.rds"))
ab5 <- readRDS(paste0(datadir, "AB5_DR.rds"))
ab6 <- readRDS(paste0(datadir, "AB6_DR.rds"))
ab7 <- readRDS(paste0(datadir, "AB7_DR.rds"))
ab8 <- readRDS(paste0(datadir, "AB8_DR.rds"))

abs <- list(ab1, ab2, ab3, ab4, ab5, ab6, ab7, ab8)

# Using 0.6 resolution for the rest of analysis
for (i in seq_along(abs)){
    ab <- abs[[i]]
    Ident(ab) <- "SCT_snn_res.0.6"
}

# Anatomical regions to annotate: Epidermis, dermis, neutrophils, fibrous capsule
# Good resource: https://ajp.amjpathol.org/article/S0002-9440%2815%2900070-X/fulltext#fig1

## AB1
ab1_neutrophils <- subset(ab1, idents = c(2))
ab1_fibrous-capsule <- subset(ab1,idents = c(1, 4, 0))
ab1_dermis <- subset(ab1, idents = c(3, 7))
ab1_epidermis <- subset(ab1, idents = c(6))

## AB2
ab2_neutrophils <- subset(ab2, idents = c(3))
ab2_fibrous-capsule <- subset(ab2, idents = c(2, 4))
ab2_dermis <- subset(ab2, idents = c(5))
ab1_epidermis <- subset(ab2, idents = c(9))

## AB3
ab3_neutrophils <- subset(ab3, idents = c(3))
ab3_fibrous-capsule <- subset(ab3, idents = c(2, 4))
ab3_dermis <- subset(ab3, idents = c(5))
ab3_epidermis <- subset(ab3, idents = c(9))

## AB4

## AB5

## AB6

## AB7

## AB8