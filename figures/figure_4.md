# Figure 4
## Load libraries
```
library("reshape2", lib.loc="~/Library/R/4.0/library")
library("ggplot2", lib.loc="~/Library/R/4.0/library")
```
## Make figure
```
# Load data
admix_win <- read.table("summary.txt", header=FALSE)
admix_win <- read.table("summary_rep.txt", header=FALSE)
counts <- read.table("window.admix.counts", header=FALSE)

# Process data
admix_merge <- (merge(admix_win, counts, by.x = c("V4","V5","V6"), by.y=c("V1","V2","V3")))
admix_7a <- melt(admix_merge,id.vars = c("V5","V6","V4","V4.y","V1"))
admix_7b <- as.data.frame(lapply(admix_7a,function(x) if(is.character(x)|is.factor(x)) gsub("V2","V4",x) else x))

# Plot tree
ggtree(tree.test, layout="rectangular",color="grey75")  %<+% key %>% collapse(node=42) +
  geom_tippoint(aes(shape=Type, fill=set_color), size=3.5, color="black") + 
  scale_shape_manual(values=c(21,24)) + 
  geom_treescale(linesize = 0.5,width=0.05, color="grey75", x = 0.5) +
  geom_tiplab(aes(label=Sample.ID.C), offset = 0.02, fontface="bold", size=2.6) +
  geom_rootedge(rootedge = 0.01, color="grey75") +
  theme(legend.position = "none") +
  scale_fill_identity()

# For each sample, plot the windowed ADMIXTURE statistics across each chromosome
A1 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_15"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
           width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))
  
A2 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_19"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                       width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A3 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_09"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A4 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_08"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A5 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_14"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A6 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_12"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A7 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_17"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))


A8 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_10"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A9 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_13"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A10 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B3_02"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A11 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_04"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                             width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A12 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_02"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A13 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_06"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A14 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_05"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A15 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_01"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A16 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_07"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A17 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="BK16_B8_11"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))

A18 <- ggplot(data=(subset(admix_7b, V4.y>0 & V1=="RT15_B3_01"))) + geom_bar(aes(x=V5, y=as.numeric(value), fill=variable, color=variable), 
                                                                              width=1, show.legend=F,stat="identity",size=0.1)  + facet_grid(.~V4,scales="free_x", space = "free_x")  +
  xlab("") + ylab("") + scale_y_continuous(expand = c(0,0), breaks=c(0.0,1.0)) + scale_x_continuous(expand = c(0,0)) + scale_fill_manual(values=c(V3="#eb3133",V4="#5ba1d9")) +scale_color_manual(values=c(V3="#eb3133",V4="#5ba1d9")) + theme(legend.position = "none") + theme_bw() +
  theme(panel.grid=element_blank(),strip.text.x = element_blank(),axis.ticks.x = element_blank(),axis.text.y = element_text(face="bold", color="black", size=8),
        axis.text.x = element_blank(),strip.text=element_text(face="bold"),axis.title.y=element_text(face="bold",size=10),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA))
```
## Merge plots
```
plot_grid(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18, ncol=1)
```
