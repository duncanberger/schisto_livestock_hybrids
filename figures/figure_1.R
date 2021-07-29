# Figure 1
```
library("reshape2", lib.loc="~/Library/R/4.0/library")
library("ggplot2", lib.loc="~/Library/R/4.0/library")
```

## Figure 1B
```
# Chromosome lenths (format: ID,start,stop)
chr <- read.csv("hy.chrom.csv", header=FALSE)

# GC content (100 kb windows)
gc <- read.table("gc_hybrid.bed", header=FALSE)

# Repeat features
reps <- read.table("rep_features.bed", header=FALSE)
LINE <- subset(reps,V5=="LINE")
LTR <- subset(reps,V5=="LTR")
PENNY <- subset(reps,V5=="PENNY")
DNA <- subset(reps,V5=="DNA")

# Coverage (100 kb windows)
cover <- read.table("coverage.txt", header=FALSE)
cover2 <- subset(cover, V4<200)

# Plot
circos.clear()
circos.par("track.height" = 0.15, start.degree = 90 , gap.degree=c(1, 1, 1, 1, 1, 1,1,16))
circos.genomicInitialize(data = chr,labels.cex = 1)
circos.track(ylim = c(0.28,  0.42), factors = gc$V1,  x=gc$V2, y=gc$V4)
circos.trackPoints(factors = gc$V1,  x=gc$V2, y=gc$V4, pch = 20, cex = 0.1,col = c('#882255','#cc6677','#ddcc77','#88ccee','#117733','#332288','#aa4499','#44aa99'))
circos.yaxis(labels.cex=0.2, side = "left", tick = T, sector.index ="1" )
circos.track(ylim = c(50,200), factors = cover2$V1,  x=cover2$V2, y=(cover2$V4))
circos.trackPoints(factors = cover2$V1,  x=cover2$V2, y=(cover2$V4), pch = 20, cex = 0.1,col = c('#882255','#cc6677','#ddcc77','#88ccee','#117733','#332288','#aa4499','#44aa99'))
circos.yaxis(labels.cex=0.2, side = "left", tick = T, sector.index ="1" )
circos.track(factors=lines$V1,x=lines$V2, y=lines$V4)
circos.trackLines(factors=LINE$V1,x=LINE$V2, y=LINE$V4, col="#CC79A7", lwd=2)
circos.trackLines(factors=LTR$V1,x=LTR$V2, y=LTR$V4, col="#E69F00", lwd=2)
circos.trackLines(factors=PENNY$V1,x=PENNY$V2, y=PENNY$V4, col="#0072B2", lwd=2)
circos.trackLines(factors=DNA$V1,x=DNA$V2, y=DNA$V4, col="#009E73", lwd=2)
circos.yaxis(labels.cex=0.2, side = "left", tick = T, sector.index ="1" )
```
## Figure 1C
```
# Read in a list of chromosome sizes
sizes <- read.table("chrom.sizes.csv", header=FALSE, sep=",",comment.char = "")

# Read in a list of promer matches
coords <- read.table("match.coords", header=FALSE)

# Reorder and filter promer hits
promer_7a <- melt(coords,id.vars = c("V5","V6","V3","V4","V9","V7","V8"))
promer_7b <- as.data.frame(lapply(promer_7a,function(x) if(is.character(x)|is.factor(x)) gsub("V1","V3",x) else x))
sizes_7a <- melt(sizes,id.vars = c("V1","V4","V5","V6","V7"))
promer_7b$x <- paste(promer_7b$V5,promer_7b$variable)

# Plot per-chromosome synteny
ggplot(data=promer_7b) +  
  geom_line(aes(x=x, y=value, group=V9, color=V5), size=0.0075) +
  geom_line(data=sizes_7a, aes(y=value, x=V5, color=as.factor(V1)),size=1.6,width = 1,lineend="round", lwd=2) +
  scale_color_manual(values=c('#882255','#cc6677','#ddcc77','#88ccee','#117733','#332288','#aa4499','#44aa99')) +
  geom_text(data=subset(sizes_7a, value!=0), aes(y=(value+6000000),x=V5, label=V7),fontface="bold", angle=90, size=3.5) +
  ylim(0,100000000) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
```
