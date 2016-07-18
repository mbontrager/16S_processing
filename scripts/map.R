# Map figure
library(ggplot2)
library(dplyr)
library(mapdata)
library(RColorBrewer)
library(gridExtra)

#setwd("manuscripts/Eaffins16S/Figs/")
sampleData <- read.csv("16S_metadata.csv", header = TRUE, 
                       na.strings = c("", "NA"))
europedata <- filter(sampleData, Continent=="Europe")
sampleData <- filter(sampleData, Continent=="North America",
                     Environment=="Copepod", 
                     !(is.na(Code)), 
                     !(Sample %in% c("VE", "V2E", "RAE", "BBE", "POE")))

redfresh <- filter(sampleData, Region=="Northeast", WaterType=="Fresh")
rf <- brewer.pal(6, "Paired")[5]
redsalt <- filter(sampleData, Region=="Northeast", WaterType=="Salt")
rs <- brewer.pal(6, "Paired")[6]
greenfresh <- filter(sampleData, Region=="South", WaterType=="Fresh")
gf <- brewer.pal(6, "Paired")[3]
greensalt <- filter(sampleData, Region=="South", WaterType=="Salt")
gs <- brewer.pal(6, "Paired")[4]

theme_nothing <- function(base_size = 12, base_family = "Helvetica")
{
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
        theme(
            line             = element_blank(),
            text             = element_blank()
        )
}

usa <- map_data("usa")
canada <- map_data("worldHires", "Canada")
NAmap <- ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = "white", 
                                        color="black") +
    geom_polygon(data = canada, aes(x=long, y = lat, group = group), fill = "white", color="black") + 
    coord_fixed(1.3) + coord_fixed(xlim = c(-100, -65),  ylim = c(25, 50), ratio = 1.3)
NAmap <- NAmap + geom_point(data=redfresh, aes(x=Longitude, y=Latitude), fill=rf, color = "black", shape=21, size=4) +
    geom_point(data=redsalt, aes(x=Longitude, y=Latitude), fill=rs, color = "black", shape=21, size=4) +
    geom_point(data=greenfresh, aes(x=Longitude, y=Latitude), fill=gf, color = "black", shape=21, size=4) +
    geom_point(data=greensalt, aes(x=Longitude, y=Latitude), fill=gs, color = "black", shape=21, size=4)
#pdf("NAmerica.pdf",width=7,height=5)
NAmap <- NAmap + theme_nothing() + theme(panel.background = element_rect(fill = "grey80")) 
#dev.off()


purpfresh <- filter(europedata, WaterType=="Fresh")
pf <- brewer.pal(10, "Paired")[9]
purpsalt <- filter(europedata, WaterType=="Salt")
ps <- brewer.pal(10, "Paired")[10]

europe <- map_data("worldHires")
EUmap <- ggplot() + geom_polygon(data = europe, aes(x=long, y = lat, group = group), 
                                 fill = "white", color="black") +
        coord_fixed(1.3) + coord_fixed(xlim = c(2, 7.6),  ylim = c(50, 54), ratio = 1.3)
EUmap <- EUmap + geom_point(data=purpfresh, aes(x=Longitude, y=Latitude), fill=pf, color = "black", shape=21, size=4) +
    geom_point(data=purpsalt, aes(x=Longitude, y=Latitude), fill=ps, color = "black", shape=21, size=4)
#pdf("Europe.pdf",width=7,height=5)
EUmap <- EUmap + theme_nothing() + theme(panel.background = element_rect(fill = "grey80")) 
#dev.off()
pdf("raw_map.pdf", width=7, height=5)
grid.arrange(NAmap, EUmap, ncol=2)
dev.off()
