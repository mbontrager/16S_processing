## Martin Bontrager
# May 13, 2015
#E. affinis 16S analysis (from mothur output)

library(phyloseq)
library(ggplot2)
library(RColorBrewer)

setwd("/home/martin/Dropbox/Projects/16S_Eaffinis/Analyses/2015-05-13/")

##Read in the .csv otu table and metadata
otuTable <- read.csv("2015-05-13_taxonomy_subtable.csv", header = TRUE, row.names = 1)
metaData <- read.csv("16S_metadata.csv", header = TRUE, row.names = 1)

# Create matrices for phyloseq
otumat <- as.matrix(subset(otuTable, select=c(8:length(colnames(otuTable)))))
otuTable[] <- lapply(otuTable, as.character)
taxmat <- as.matrix(subset(otuTable, select=c(1:7)))

## Declare phyloseq objects
OTU <- otu_table(otumat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
metaData <- sample_data(metaData)

physeq <- phyloseq(OTU, TAX, metaData)

## GGplot themeing
theme_set(theme_bw())
pal = "Set1"
scale_colour_discrete <- function(palname = pal, ...) {
        scale_colour_brewer(palette = palname, ...)
}
scale_fill_discrete <- function(palname = pal, ...) {
        scale_fill_brewer(palette = palname, ...)
}

#Prune taxa not present in any sample (if they exist)
GP <- prune_taxa(taxa_sums(physeq) > 0, physeq)

#Generate a large richness plot(many different estimators)
plot_richness(GP)

# Plot Richness only with chao1 and shannon:
plot_richness(GP,measures = c("Chao1", "Shannon"))

# Plot richness Salt water vs. Freshwater
plot_richness(GP, x = "WaterType", measures = c("Chao1", "Shannon"))

# Plot richness copepod vs. water:
colorCount <- nsamples(GP)
getPalette <- colorRampPalette(brewer.pal(8, "Set1"))
p <- plot_richness(GP, x = "Environment", color = "Location", 
              measures = c("Chao1", "Shannon")) + 
        scale_colour_manual(values = getPalette(colorCount))
p <- p + geom_point(size = 3, alpha = 0.8)
p

#Merge samples by region and plot richness again:
GP_merge1 <- merge_samples(GP,"Region")
colorCount <- nsamples(GP_merge1)
plot_richness(GP_merge1, x = "Environment", color = "Region", 
              measures = c("Chao1", "Shannon")) + 
        scale_colour_manual(values = getPalette(colorCount))

# Plot richness copepod vs. water:
colorCount <- nsamples(GP)
getPalette <- colorRampPalette(brewer.pal(8, "Set1"))
p <- plot_richness(GP, x = "Environment", color = "Region", 
              measures = c("Chao1", "Shannon"))
             
p <- p + geom_point(size = 5, alpha = 0.7)

##Bar Plots of taxa abundance
##Just the proteobacteria
GP_sub1 <- subset_taxa(GP, Phylum == "Proteobacteria")
plot_bar(GP_sub1)
