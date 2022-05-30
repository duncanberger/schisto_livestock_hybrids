## Figure 2
### Load libraries
```
library("ggplot2")
library("ggrepel")
library("cowplot")
library("ape", lib.loc="~/Library/R/4.0/library")
library("phytools", lib.loc="~/Library/R/4.0/library")
library("ggtree", lib.loc="~/Library/R/4.0/library")
```
### Autosomal PCA
```
# Load metadata
key <- read.table("metadata.csv", header=TRUE, sep=",",comment.char = "")

# Load plink PCA files
df<-read.delim("prunedData.eigenvec", header=TRUE, sep="\t")
df_2 <- (merge(key, df, all=TRUE, by.y = "IID", by.x='sample'))
eigens<-read.delim("prunedData.eigenval",header=F)
sum_eigs<-sum(eigens$V1)
sum_eigs<-lapply(eigens$V1,function(x){
  rt<-(x/sum_eigs)*100
  rt<-round(rt)
  return(rt)
})

# Plot PCA
pca1 <-  ggplot(df_2, aes((PC1),(PC2))) +
  geom_point(aes(fill=set_color, shape=Type),size = 2.5, alpha=0.7) +
  scale_shape_manual(values=c(22,21,24)) +
  xlab(paste0("PC1 (",sum_eigs[[1]],"%)")) + ylab(paste0("PC2 (",sum_eigs[[2]],"%)")) +
  theme_bw() + 
  scale_x_continuous(limits=c(-0.3,0.3), expand=c(0,0), breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3)) +
  scale_y_continuous(limits=c(-0.3,0.7), expand=c(0,0), breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
  scale_fill_identity() +
  theme(panel.grid=element_blank(), 
        axis.title=element_text(face="bold", size=7),
        legend.text = element_text(face="bold"), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=5, color="black")) 

pca3 <- ggplot(df_2, aes((PC3),(PC4))) +
  geom_point(aes(fill=set_color, shape=Type),size = 2.5, alpha=0.7) +
  scale_shape_manual(values=c(22,21,24)) +
  xlab(paste0("PC3 (",sum_eigs[[3]],"%)")) + ylab(paste0("PC4 (",sum_eigs[[4]],"%)")) +
  theme_bw() + 
  scale_fill_identity() +
  scale_x_continuous(limits=c(-0.3,0.4), expand=c(0,0), breaks=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4)) +
  scale_y_continuous(limits=c(-0.5,0.7), expand=c(0,0), breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)) +
  theme(panel.grid=element_blank(), 
        axis.title=element_text(face="bold", size=7),
        legend.text = element_text(face="bold"), 
        legend.position = "none", 
        axis.text = element_text(face = "bold", size=5)) 
```
### Phylogeny
```
# Load ML tree
treex.test <- midpoint(read.newick("ML_wgs.nwk"))
auto2 <- (root(treex.test, node=37))

# Plot tree
ggtree(tree.test, layout="ape",color="grey75")  %<+% key+
  geom_tippoint(aes(shape=Type, fill=set_color), size=2, color="black") + 
  scale_shape_manual(values=c(21,24)) +
  geom_treescale(linesize = 0.5,width=0.05, color="grey75") +
  theme(legend.position = "none") + scale_colour_identity() + 
  scale_fill_identity() 
```
### ADMIXTURE
```
# Load data
admix <- read.table("admixture_all.txt", sep="\t", header=FALSE)
admix_summary <- (melt(admix,id.vars = c("V1","V5")))
admix_summary_merged <- (merge(key, admix_summary, by.y = "V1", by.x='sample'))

# Plot
admix <- ggplot(data=(subset(admix_summary_merged)) +
  geom_bar(aes(x=V5, y=as.numeric(value), fill=variable), 
           width=1, show.legend=F,stat="identity", color="black",size=0.1)  + 
  facet_grid(.~Group,scales="free_x", space = "free_x") + 
  xlab("") + ylab("Admixture proportion") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c(V3="#eb3133",V4="#984ea3",V2="#5ba1d9")) +theme(legend.position = "none") +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face="bold", color="black", size=5),
        axis.text.x = element_blank(),
        strip.text=element_text(face="bold"),
        axis.title.y=element_text(face="bold",size=7),
        panel.background = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color="black",fill=NA))
```
#### Merge plots
```
plot1 <- plot_grid(pca1,pca3, ncol=1, align="v", label_size = 10, labels=c("A","B") )
plot2 <- plot_grid(plot1,tree,tree, ncol=3, labels=c("","C",""), rel_widths = c(0.8,1,0.1), label_size = 10)
plot_grid(plot2,admix, nrow=2, labels=c("","D"), rel_heights = c(1,0.35), label_size = 10)
```
