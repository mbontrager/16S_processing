## Plot an ordered bar chart (easier for grouping samples non-alphabetically)
## Use the ordered plot function to keep phyla in a consistent position
# Only works on unfaceted bar plots.

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