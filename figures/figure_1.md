# Figure 1
```
library("reshape2", lib.loc="~/Library/R/4.0/library")
library("ggplot2", lib.loc="~/Library/R/4.0/library")
library("circlize", lib.loc="~/Library/R/4.0/library")
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
