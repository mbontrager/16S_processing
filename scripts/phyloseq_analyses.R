## 16S analysis with Phyloseq and ggplot2. Installation instructions are online
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
theme_set(theme_bw())

myData = import_biom("D:/microbiome5-1-14/otu.biom", parseFunction = parse_taxonomy_greengenes)
plot_bar(myData, fill = "Class")
