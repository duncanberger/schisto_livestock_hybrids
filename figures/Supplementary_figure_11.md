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
tree.test <- midpoint(read.tree("mito_phylo.nwk"))

# Plot tree
tree.test2 <- root(tree.test, node = 40)
mito_phylo_1 <- ggtree(tree.test, layout = "rectangular",color="grey75")  %<+% key +
  geom_tippoint(aes(shape=Type, fill=set_color), size=2.5, color="black") + 
  scale_shape_manual(values=c(21,24)) +
  geom_text2(aes(subset=!isTip, label=label), size=1.95, offset = 0.025) +
  geom_tiplab(aes(label=Sample.ID.C), fontface="bold", size=2, offset = 0.025) +
  geom_treescale(linesize = 0.5,width=0.05, color="grey75", x=0.75) +
  theme(legend.position = "none") + scale_colour_identity() + 
  scale_fill_identity() 

mito_phylo_2 <- ggtree(tree.test, layout = "ape",color="grey75")  %<+% key +
  geom_tippoint(aes(shape=Type, fill=set_color), size=2.5, color="black") + 
  scale_shape_manual(values=c(21,24)) +
  geom_treescale(linesize = 0.5,width=0.05, color="grey75") +
  theme(legend.position = "none") + scale_colour_identity() + 
  scale_fill_identity() 
```
