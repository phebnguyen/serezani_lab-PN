#Import packages
library(ggplot2)

##### Configuration #####
setwd("/Users/nguyenpu/Documents/GitHub/serezani_lab/bulk_rna_seq/pten_ko/")
datadir <- "outputdir/"
outputdir <- "outputdir/volcano_plots"

if (dir.exists(outputdir) != TRUE) {
    dir.create(outputdir)
}

#DESeq2 Data
uninfected_deseq2 <- read.csv(paste0(datadir, "DESeq2_CTLUninfected_vs_KOUninfected.csv"))
infected_deseq2 <- read.csv(paste0(datadir, "DESeq2_CTLInfected_vs_KOInfected.csv"))

##### Plotting #####
#Write basic volcano plot function
volcano_plot <- function(deseq2_data, file_name, plot_title, outputdir){
    color_map <- c("UP" = "red", "DOWN" = "blue", "NO" = "grey")
    deseq2_data$padj[deseq2_data$padj == 0] <- 1e-300
    upper_limit <- max(-log10(deseq2_data$padj), na.rm = TRUE) * 1.05  # Add 10% space

    p <- ggplot(data = deseq2_data, aes(x = log2FoldChange, y = -log10(padj), color = DEGenes)) + 
        geom_point() + 
        scale_color_manual(values = color_map) +
        theme_classic() + 
        xlab(expression(log[2] ~ "Fold Change")) +
        ylab(expression(-log[10] ~ "Adjusted P-value")) +
        ggtitle(plot_title) + 
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_text(size = 14),
              axis.title.y = element_text(size = 14),
              plot.title = element_text(size = 16),
              legend.position = "none") +
        geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dotted") +
        geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted") +
        coord_cartesian(ylim = c(0, upper_limit))

    # Save the plot
    ggsave(filename = paste0(outputdir, "/", file_name), plot = p, 
           width = 10, height = 10, units = "in", bg = "white")
}

#Create volcano plots for uninfected and infected data
volcano_plot(uninfected_deseq2, "uninfected_volcano_plot.png", "Uninfected KO vs. WT", outputdir)
volcano_plot(infected_deseq2, "infected_volcano_plot.png", "Infected KO vs. WT", outputdir)
