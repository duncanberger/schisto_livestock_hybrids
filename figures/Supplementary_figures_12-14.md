# Supplementary figures 12-14
## Load libraries
```
library("ggplot2")
```
## Read in data
```
sw <- read.table("median.sliding.windows.txt", header=FALSE)
sz <- read.csv("sizes_windows.csv", header=FALSE)
```
## Plot for each sample
```
X1 <- ggplot(data=subset(sw, V4=="BK16_B8_15")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X2 <- ggplot(data=subset(sw, V4=="BK16_B8_19")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X3 <- ggplot(data=subset(sw, V4=="BK16_B8_09")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X4 <- ggplot(data=subset(sw, V4=="BK16_B8_08")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X5 <- ggplot(data=subset(sw, V4=="BK16_B8_14")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X6 <- ggplot(data=subset(sw, V4=="BK16_B8_12")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X7 <- ggplot(data=subset(sw, V4=="BK16_B8_17")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X8 <- ggplot(data=subset(sw, V4=="BK16_B8_10")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X9 <- ggplot(data=subset(sw, V4=="BK16_B8_13")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X10 <- ggplot(data=subset(sw, V4=="BK16_B3_02")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X11 <- ggplot(data=subset(sw, V4=="BK16_B8_04")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X12 <- ggplot(data=subset(sw, V4=="BK16_B8_02")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X13 <- ggplot(data=subset(sw, V4=="BK16_B8_06")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X14 <- ggplot(data=subset(sw, V4=="BK16_B8_05")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X15 <- ggplot(data=subset(sw, V4=="BK16_B8_01")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X16 <- ggplot(data=subset(sw, V4=="BK16_B8_07")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X17 <- ggplot(data=subset(sw, V4=="BK16_B8_11")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")

X18 <- ggplot(data=subset(sw, V4=="RT15_B3_01")) + 
  geom_point(aes(x=V2/1000000,y=V5), size=0.1)  + 
  geom_line(aes(x=V2/1000000,y=V5), size=0.1)  + 
  facet_grid(.~V1,scales="free_x", space = "free_x") +
  geom_rect(data=sz, inherit.aes=FALSE, aes(xmin=V2/1000000, xmax=V3/1000000, ymin=V4,ymax=V5, fill=V6), alpha=0.2) +
  scale_fill_manual(values=c("#5ba1d9","#eb3133")) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous( limits=c(0,1))+
  ylab("Admixture proportion") + xlab("Position on chromosome (Mb)") +
  theme(panel.grid=element_blank(), axis.text.x =element_text(face="bold", color="black", size=5),
        axis.title.x =element_text(face="bold", color="black", size=6),
        axis.text.y = element_text(face="bold", color="black", size=5),
        strip.text = element_text(face="bold", color="black", size=5),
        axis.title.y=element_text(face="bold",size=6),panel.background = element_blank(),
        strip.background = element_blank(),panel.border = element_rect(color="black",fill=NA), legend.position = "none")
```
# Merge figures
```
XX01<- plot_grid(X10,X15,X12,X11,X14,X13, ncol=1, 
                 labels=c('BK16_B3_02','BK16_B8_01','BK16_B8_02',"BK16_B8_04",
                          "BK16_B8_05","BK16_B8_06"), 
                 label_size = 8)

XX02 <- plot_grid(X16,X4,X3,X8,X17,X6, ncol=1, 
                 labels=c("BK16_B8_07","BK16_B8_08","BK16_B8_09","BK16_B8_10","BK16_B8_11","BK16_B8_12"), 
                 label_size = 8)

XX03 <- plot_grid(X9,X5,X1,X7,X2,X18, ncol=1, 
                  labels=c("BK16_B8_13","BK16_B8_14","BK16_B8_15","BK16_B8_17","BK16_B8_19","RT15_B3_01"), 
                  label_size = 8)
```
