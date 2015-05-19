## Martin Bontrager
# May 13, 2015
#E. affinis 16S analysis (from mothur output)

library(phyloseq)
library(ggplot2)
library(RColorBrewer)

setwd("D:/Users/Martin/Dropbox/Projects/16S_Eaffinis/Analyses/2015-05-13/")

##Read in the .csv otu table and metadata
otuTable <- read.csv("2015-05-13_taxonomy_subtable.csv", header = TRUE, row.names = 1)
metaData <- read.csv("16S_metadata.csv", header = TRUE, row.names = 1)
ordered_samples <- as.vector((read.table("samples.txt"))[, 1])

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
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

#Prune taxa not present in any sample (if they exist)
#And remove unnecessary objects from memory
GP <- prune_taxa(taxa_sums(physeq) > 0, physeq)
rm(otuTable, otumat, taxmat, physeq, TAX, OTU)

#Generate a large richness plot(many different estimators)
plot_richness(GP)

# Plot Richness only with chao1 and shannon:
plot_richness(GP,measures = c("Chao1", "Shannon"))

# Plot richness Salt water vs. Freshwater
plot_richness(GP, x = "WaterType", measures = c("Chao1", "Shannon"))

# Plot richness copepod vs. water:
colorCount <- nsamples(GP)
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
getPalette <- colorRampPalette(brewer.pal(9, "Set2"))
p <- plot_richness(GP, x = "Environment", color = "Region", 
              measures = c("Chao1", "Shannon"))
             
p <- p + geom_point(size = 5, alpha = 0.7)

##Bar Plots of taxa abundance
## Overall
##Find out the number of Phyla, for example:
a <- as.vector(tax_table(GP)[, "Phylum"])
b <- unique(a)
rm(a)
colorCount <- length(b)

# Merge (glom) the physeq object down to phylum level for bar plot gen:
##Must use "scale_fill_manual" for fill objects instead of `scale_color_manual`
##which is used for points/lines

GP_glom1 <- tax_glom(GP, "Phylum")
p <- plot_bar(GP_glom1, fill="Phylum") +
        scale_fill_manual(values = getPalette(colorCount)) +
        guides(fill=guide_legend(ncol=2))

## Use the ordered plot function to keep phyla in a consistent position
# leg_size scales the size of the legend (use 1.0 for large)
# scale_x_disctrete uses an ordered list of samples that I read in
# to organize the x axis instead of alphabetical order.
p1 <- plot_ordered_bar(GP_glom1, fill="Phylum", leg_size = 0.8) +
        scale_fill_manual(values = getPalette(colorCount)) + 
        guides(fill=guide_legend(ncol=2)) +
        scale_x_discrete(limits = ordered_samples)

##Just the Flavobacterium
GP_sub1 <- subset_taxa(GP, Genus == "Flavobacterium")
theme_set(theme_bw() +
                  theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()))
plot_bar(GP_sub1, fill = "Species", title = "Flavobacterium") + 
        scale_x_discrete(limits = ordered_samples) + 
        theme(axis.title.x = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              plot.title = element_text(size = 24, face = "bold"))

##This function will plot an ordered bar chart, consistent across samples:
plot_ordered_bar<-function (physeq, x = "Sample", 
                            y = "Abundance", 
                            fill = NULL, 
                            leg_size = 0.5,
                            title = NULL, facet_grid = NULL) 
{
        require(ggplot2)
        require(phyloseq)
        require(plyr)
        require(grid)
        bb <- psmelt(physeq)
        .e <- environment()
        p = ggplot(bb, aes_string(x = x, y = y, 
                                  fill = fill), environment = .e)
        
        
        p = p + geom_bar(aes(order=desc(bb[,fill])),
                         stat = "identity", 
                         position = "stack", 
                         color = "black") 
        
        p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
        
        p = p + guides(fill = guide_legend(override.aes = list(colour = NULL))) 
        p = p + theme(legend.key = element_rect(colour = "black")) 
        
        p = p + theme(legend.key.size = unit(leg_size, "cm"))
        
        if (!is.null(facet_grid)) {
                
                p <- p + facet_grid(facet_grid)
        }
        
        if (!is.null(title)) {
                p <- p + ggtitle(title)
        }
        return(p)
}
