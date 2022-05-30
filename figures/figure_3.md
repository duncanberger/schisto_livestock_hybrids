# Figure 3
## Load libraries
```
library("ggplot2")
library("ape")
library("ggtree")
library(see)
library("ggstance")
library("reshape2")
library("dplyr", lib.loc="~/Library/R/4.0/library")
library("treeio", lib.loc="~/Library/R/4.0/library")
```
### Plot Figure
```
# Read in metadata
key <- read.table("metadata.csv", header=TRUE, sep=",",comment.char = "")

# Read in data
hets <- read.table("merged.hets.all.fixed2txt", header=FALSE)
admix <- read.table("all.admix.final", sep='\t')
f3 <- read.table("f3.refixed.txt", header=FALSE)
tree.test <- read.newick("auto_phylo.nwk")
D <- read.table("d.fixed.txt", sep="\t")

# Process data
hets_subset <- hets %>% group_by(V1) %>% sample_n(size=5000)
f3_subset <- f3 %>% group_by(V1,V2) %>% sample_n(size=5000)
f3_subset_2 <- subset(f3_subset, V3<0.25)
D_subset <- D %>% group_by(V1,V2) %>% sample_n(size=5000)
admix_2 <- (melt(admix,id.vars = c("V1","V5")))

# Plot tree
b <- ggtree(tree.test, layout="rectangular",color="grey75")  %<+% key %>% collapse(node=58) +
  geom_tippoint(aes(shape=Type, fill=set_color), size=1.5, color="black") + 
  scale_shape_manual(values=c(21,24)) + 
  geom_treescale(linesize = 0.25,width=0.05, color="grey75") +
  geom_tiplab(aes(label=Sample.ID.C), offset = 0.04, fontface="bold", size=1.5) +
  geom_rootedge(rootedge = 0.01, color="grey75") +
  theme(legend.position = "none") +
  scale_fill_identity()

# Plot ADMIXTURE
c <- facet_plot(b, panel="adm",
                data=admix_2,
                geom=geom_colh,
                mapping=aes(x=as.numeric(value), fill=set_color) , size=0.1) + 
  theme_tree2() + 
  theme(legend.position = "none")

# Plot heterozygous SNPs
d <- facet_plot(c, panel="HET",
                data=hets_subset,
                geom=geom_jitter2,
                mapping=aes(x=((V6/V5)*100), color="grey50"),fill="grey50", size=0.15, alpha=0.4) +
  scale_color_identity() 

e <- facet_plot(d, panel="HET",
           data=hets_subset,
           geom=geom_boxploth,
           mapping=aes(x=((V6/V5)*100), group=V7),alpha=0,notch=TRUE, outlier.shape=NA) 

# Plot f3 statistics
f <- facet_plot(e, panel="f3",
                data=f3_subset_2,
                geom=geom_jitter2,
                mapping=aes(x=(V3), group=V2, color="grey50"),fill="grey50", size=0.15, alpha=0.4)

j <- facet_plot(f, panel="f3",
                data=(f3_subset_2),
                geom=geom_boxploth,
                mapping=aes(x=(V3), group=V2),alpha=0,notch=TRUE, outlier.shape=NA)

# Resize tree
j<- i + xlim_tree(0.5) + 
  theme(axis.text.x = element_text(color="black", size=5))

# Save image to file
ggsave(filename = "fig3.svg",j, units = c("cm"), width = 16.9, height=8)
```
