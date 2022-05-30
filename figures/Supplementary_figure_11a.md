# Supplementary figure 11
## Load libraries
```
library("ggplot2")
library("cowplot")
library("ape", lib.loc="~/Library/R/4.0/library")
library("phytools", lib.loc="~/Library/R/4.0/library")
library("ggtree", lib.loc="~/Library/R/4.0/library")
```
## Mitochondrial phylogeny
```
# Read tree
tree.test <- read.newick("mito_ML_tree.nwk")

# Plot tree
mito_tree <- ggtree(tree.test, layout = "ape", color="grey75")  %<+% key +
  geom_tippoint(aes(shape=Type, fill=set_color), size=2.5, color="black") + 
  scale_shape_manual(values=c(21,24)) +
  geom_treescale(linesize = 0.5,width=0.05, color="grey75") +
  theme(legend.position = "none") + scale_colour_identity() + 
  scale_fill_identity() 
```
