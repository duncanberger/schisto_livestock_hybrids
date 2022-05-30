# Supplementary figures 1-4
## Load libraries
```
library("ggplot2")
library("reshape2")
```
## Load data and formatting
```
coverage <- read.table("all.recov.txt", header=FALSE)
selection_colors <- rep(c("grey75", "grey40"))
selection_theme <- theme(
  legend.position="none",
  panel.grid = element_blank(),
  axis.text.x=element_text(face="bold", color="black"),
  axis.text.y=element_text(face="bold", color="black", size=5),
  axis.title.y = element_text(face="bold", color="black", size=5),
  axis.title.x = element_blank(),
  axis.ticks.x=element_blank())
```
## Make plots (for each sample, using internal sample names)
```
cov_hybrid_schisto7118176b <- subset(coverage, V7=="hybrid_schisto7118176b")
cov_hybrid_schisto7118176b<- cov_hybrid_schisto7118176b[order(cov_hybrid_schisto7118176b$V1, cov_hybrid_schisto7118176b$V2),]
cov_hybrid_schisto7118176b$ID <- seq.int(nrow(cov_hybrid_schisto7118176b))

axisdf = cov_hybrid_schisto7118176b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B3_01.out <- ggplot(subset(cov_hybrid_schisto7118176b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118177b <- subset(coverage, V7=="hybrid_schisto7118177b")
cov_hybrid_schisto7118177b<- cov_hybrid_schisto7118177b[order(cov_hybrid_schisto7118177b$V1, cov_hybrid_schisto7118177b$V2),]
cov_hybrid_schisto7118177b$ID <- seq.int(nrow(cov_hybrid_schisto7118177b))

axisdf = cov_hybrid_schisto7118177b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B3_02.out <- ggplot(subset(cov_hybrid_schisto7118177b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118178b <- subset(coverage, V7=="hybrid_schisto7118178b")
cov_hybrid_schisto7118178b<- cov_hybrid_schisto7118178b[order(cov_hybrid_schisto7118178b$V1, cov_hybrid_schisto7118178b$V2),]
cov_hybrid_schisto7118178b$ID <- seq.int(nrow(cov_hybrid_schisto7118178b))

axisdf = cov_hybrid_schisto7118178b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_01.out <- ggplot(subset(cov_hybrid_schisto7118178b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118179b <- subset(coverage, V7=="hybrid_schisto7118179b")
cov_hybrid_schisto7118179b<- cov_hybrid_schisto7118179b[order(cov_hybrid_schisto7118179b$V1, cov_hybrid_schisto7118179b$V2),]
cov_hybrid_schisto7118179b$ID <- seq.int(nrow(cov_hybrid_schisto7118179b))

axisdf = cov_hybrid_schisto7118179b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_02.out <- ggplot(subset(cov_hybrid_schisto7118179b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118180b <- subset(coverage, V7=="hybrid_schisto7118180b")
cov_hybrid_schisto7118180b<- cov_hybrid_schisto7118180b[order(cov_hybrid_schisto7118180b$V1, cov_hybrid_schisto7118180b$V2),]
cov_hybrid_schisto7118180b$ID <- seq.int(nrow(cov_hybrid_schisto7118180b))

axisdf = cov_hybrid_schisto7118180b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_03.out <- ggplot(subset(cov_hybrid_schisto7118180b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118181b <- subset(coverage, V7=="hybrid_schisto7118181b")
cov_hybrid_schisto7118181b<- cov_hybrid_schisto7118181b[order(cov_hybrid_schisto7118181b$V1, cov_hybrid_schisto7118181b$V2),]
cov_hybrid_schisto7118181b$ID <- seq.int(nrow(cov_hybrid_schisto7118181b))

axisdf = cov_hybrid_schisto7118181b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_04.out <- ggplot(subset(cov_hybrid_schisto7118181b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118182b <- subset(coverage, V7=="hybrid_schisto7118182b")
cov_hybrid_schisto7118182b<- cov_hybrid_schisto7118182b[order(cov_hybrid_schisto7118182b$V1, cov_hybrid_schisto7118182b$V2),]
cov_hybrid_schisto7118182b$ID <- seq.int(nrow(cov_hybrid_schisto7118182b))

axisdf = cov_hybrid_schisto7118182b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_05.out <- ggplot(subset(cov_hybrid_schisto7118182b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7118183b <- subset(coverage, V7=="hybrid_schisto7118183b")
cov_hybrid_schisto7118183b<- cov_hybrid_schisto7118183b[order(cov_hybrid_schisto7118183b$V1, cov_hybrid_schisto7118183b$V2),]
cov_hybrid_schisto7118183b$ID <- seq.int(nrow(cov_hybrid_schisto7118183b))

axisdf = cov_hybrid_schisto7118183b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_06.out <- ggplot(subset(cov_hybrid_schisto7118183b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7148988b <- subset(coverage, V7=="hybrid_schisto7148988b")
cov_hybrid_schisto7148988b<- cov_hybrid_schisto7148988b[order(cov_hybrid_schisto7148988b$V1, cov_hybrid_schisto7148988b$V2),]
cov_hybrid_schisto7148988b$ID <- seq.int(nrow(cov_hybrid_schisto7148988b))

axisdf = cov_hybrid_schisto7148988b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
RT15_B6_01.out <- ggplot(subset(cov_hybrid_schisto7148988b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714151b <- subset(coverage, V7=="hybrid_schisto7714151b")
cov_hybrid_schisto7714151b<- cov_hybrid_schisto7714151b[order(cov_hybrid_schisto7714151b$V1, cov_hybrid_schisto7714151b$V2),]
cov_hybrid_schisto7714151b$ID <- seq.int(nrow(cov_hybrid_schisto7714151b))

axisdf = cov_hybrid_schisto7714151b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_07.out <- ggplot(subset(cov_hybrid_schisto7714151b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714152b <- subset(coverage, V7=="hybrid_schisto7714152b")
cov_hybrid_schisto7714152b<- cov_hybrid_schisto7714152b[order(cov_hybrid_schisto7714152b$V1, cov_hybrid_schisto7714152b$V2),]
cov_hybrid_schisto7714152b$ID <- seq.int(nrow(cov_hybrid_schisto7714152b))

axisdf = cov_hybrid_schisto7714152b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_08.out <- ggplot(subset(cov_hybrid_schisto7714152b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714153b <- subset(coverage, V7=="hybrid_schisto7714153b")
cov_hybrid_schisto7714153b<- cov_hybrid_schisto7714153b[order(cov_hybrid_schisto7714153b$V1, cov_hybrid_schisto7714153b$V2),]
cov_hybrid_schisto7714153b$ID <- seq.int(nrow(cov_hybrid_schisto7714153b))

axisdf = cov_hybrid_schisto7714153b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_09.out <- ggplot(subset(cov_hybrid_schisto7714153b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714154b <- subset(coverage, V7=="hybrid_schisto7714154b")
cov_hybrid_schisto7714154b<- cov_hybrid_schisto7714154b[order(cov_hybrid_schisto7714154b$V1, cov_hybrid_schisto7714154b$V2),]
cov_hybrid_schisto7714154b$ID <- seq.int(nrow(cov_hybrid_schisto7714154b))

axisdf = cov_hybrid_schisto7714154b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_10.out <- ggplot(subset(cov_hybrid_schisto7714154b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,125)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714155b <- subset(coverage, V7=="hybrid_schisto7714155b")
cov_hybrid_schisto7714155b<- cov_hybrid_schisto7714155b[order(cov_hybrid_schisto7714155b$V1, cov_hybrid_schisto7714155b$V2),]
cov_hybrid_schisto7714155b$ID <- seq.int(nrow(cov_hybrid_schisto7714155b))

axisdf = cov_hybrid_schisto7714155b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_11.out <- ggplot(subset(cov_hybrid_schisto7714155b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714156b <- subset(coverage, V7=="hybrid_schisto7714156b")
cov_hybrid_schisto7714156b<- cov_hybrid_schisto7714156b[order(cov_hybrid_schisto7714156b$V1, cov_hybrid_schisto7714156b$V2),]
cov_hybrid_schisto7714156b$ID <- seq.int(nrow(cov_hybrid_schisto7714156b))

axisdf = cov_hybrid_schisto7714156b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_12.out <- ggplot(subset(cov_hybrid_schisto7714156b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714157b <- subset(coverage, V7=="hybrid_schisto7714157b")
cov_hybrid_schisto7714157b<- cov_hybrid_schisto7714157b[order(cov_hybrid_schisto7714157b$V1, cov_hybrid_schisto7714157b$V2),]
cov_hybrid_schisto7714157b$ID <- seq.int(nrow(cov_hybrid_schisto7714157b))

axisdf = cov_hybrid_schisto7714157b %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_13.out <- ggplot(subset(cov_hybrid_schisto7714157b), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714158 <- subset(coverage, V7=="hybrid_schisto7714158")
cov_hybrid_schisto7714158<- cov_hybrid_schisto7714158[order(cov_hybrid_schisto7714158$V1, cov_hybrid_schisto7714158$V2),]
cov_hybrid_schisto7714158$ID <- seq.int(nrow(cov_hybrid_schisto7714158))

axisdf = cov_hybrid_schisto7714158 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_14.out <- ggplot(subset(cov_hybrid_schisto7714158), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714159 <- subset(coverage, V7=="hybrid_schisto7714159")
cov_hybrid_schisto7714159<- cov_hybrid_schisto7714159[order(cov_hybrid_schisto7714159$V1, cov_hybrid_schisto7714159$V2),]
cov_hybrid_schisto7714159$ID <- seq.int(nrow(cov_hybrid_schisto7714159))

axisdf = cov_hybrid_schisto7714159 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_15.out <- ggplot(subset(cov_hybrid_schisto7714159), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714161 <- subset(coverage, V7=="hybrid_schisto7714161")
cov_hybrid_schisto7714161<- cov_hybrid_schisto7714161[order(cov_hybrid_schisto7714161$V1, cov_hybrid_schisto7714161$V2),]
cov_hybrid_schisto7714161$ID <- seq.int(nrow(cov_hybrid_schisto7714161))

axisdf = cov_hybrid_schisto7714161 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_17.out <- ggplot(subset(cov_hybrid_schisto7714161), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7714163 <- subset(coverage, V7=="hybrid_schisto7714163")
cov_hybrid_schisto7714163<- cov_hybrid_schisto7714163[order(cov_hybrid_schisto7714163$V1, cov_hybrid_schisto7714163$V2),]
cov_hybrid_schisto7714163$ID <- seq.int(nrow(cov_hybrid_schisto7714163))

axisdf = cov_hybrid_schisto7714163 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
BK16_B8_19.out <- ggplot(subset(cov_hybrid_schisto7714163), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,125)) +
  theme_bw() +
  selection_theme
cov_hybrid_schisto7717246 <- subset(coverage, V7=="hybrid_schisto7717246")
cov_hybrid_schisto7717246<- cov_hybrid_schisto7717246[order(cov_hybrid_schisto7717246$V1, cov_hybrid_schisto7717246$V2),]
cov_hybrid_schisto7717246$ID <- seq.int(nrow(cov_hybrid_schisto7717246))

axisdf = cov_hybrid_schisto7717246 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
RT15_B3_01.out <- ggplot(subset(cov_hybrid_schisto7717246), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SBOV_01 <- subset(coverage, V7=="SBOV_1")
cov_SBOV_01<- cov_SBOV_01[order(cov_SBOV_01$V1, cov_SBOV_01$V2),]
cov_SBOV_01$ID <- seq.int(nrow(cov_SBOV_01))

axisdf = cov_SBOV_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SBOV_01.out <- ggplot(subset(cov_SBOV_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SCUR_01 <- subset(coverage, V7=="SCUR_01")
cov_SCUR_01<- cov_SCUR_01[order(cov_SCUR_01$V1, cov_SCUR_01$V2),]
cov_SCUR_01$ID <- seq.int(nrow(cov_SCUR_01))

axisdf = cov_SCUR_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SCUR_01.out <- ggplot(subset(cov_SCUR_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SGUI_01 <- subset(coverage, V7=="SGUI_01")
cov_SGUI_01<- cov_SGUI_01[order(cov_SGUI_01$V1, cov_SGUI_01$V2),]
cov_SGUI_01$ID <- seq.int(nrow(cov_SGUI_01))

axisdf = cov_SGUI_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SGUI_01.out <- ggplot(subset(cov_SGUI_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SHAE_01 <- subset(coverage, V7=="SHAE_01b")
cov_SHAE_01<- cov_SHAE_01[order(cov_SHAE_01$V1, cov_SHAE_01$V2),]
cov_SHAE_01$ID <- seq.int(nrow(cov_SHAE_01))

axisdf = cov_SHAE_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SHAE_01.out <- ggplot(subset(cov_SHAE_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SINT_01 <- subset(coverage, V7=="SINT_01")
cov_SINT_01<- cov_SINT_01[order(cov_SINT_01$V1, cov_SINT_01$V2),]
cov_SINT_01$ID <- seq.int(nrow(cov_SINT_01))

axisdf = cov_SINT_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SINT_01.out <- ggplot(subset(cov_SINT_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SMAG_01 <- subset(coverage, V7=="SMAG_01")
cov_SMAG_01<- cov_SMAG_01[order(cov_SMAG_01$V1, cov_SMAG_01$V2),]
cov_SMAG_01$ID <- seq.int(nrow(cov_SMAG_01))

axisdf = cov_SMAG_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SMAG_01.out <- ggplot(subset(cov_SMAG_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SMAT_01 <- subset(coverage, V7=="SMAT_01")
cov_SMAT_01<- cov_SMAT_01[order(cov_SMAT_01$V1, cov_SMAT_01$V2),]
cov_SMAT_01$ID <- seq.int(nrow(cov_SMAT_01))

axisdf = cov_SMAT_01 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SMAT_01.out <- ggplot(subset(cov_SMAT_01), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SBOV_02 <- subset(coverage, V7=="SBOV_02")
cov_SBOV_02<- cov_SBOV_02[order(cov_SBOV_02$V1, cov_SBOV_02$V2),]
cov_SBOV_02$ID <- seq.int(nrow(cov_SBOV_02))

axisdf = cov_SBOV_02 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SBOV_02.out <- ggplot(subset(cov_SBOV_02), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SHAE_02 <- subset(coverage, V7=="SHAE_02")
cov_SHAE_02<- cov_SHAE_02[order(cov_SHAE_02$V1, cov_SHAE_02$V2),]
cov_SHAE_02$ID <- seq.int(nrow(cov_SHAE_02))

axisdf = cov_SHAE_02 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SHAE_02.out <- ggplot(subset(cov_SHAE_02), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SHAE_04 <- subset(coverage, V7=="SHAE_04")
cov_SHAE_04<- cov_SHAE_04[order(cov_SHAE_04$V1, cov_SHAE_04$V2),]
cov_SHAE_04$ID <- seq.int(nrow(cov_SHAE_04))

axisdf = cov_SHAE_04 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SHAE_04.out <- ggplot(subset(cov_SHAE_04), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SHAE_05 <- subset(coverage, V7=="SHAE_05")
cov_SHAE_05<- cov_SHAE_05[order(cov_SHAE_05$V1, cov_SHAE_05$V2),]
cov_SHAE_05$ID <- seq.int(nrow(cov_SHAE_05))

axisdf = cov_SHAE_05 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SHAE_05.out <- ggplot(subset(cov_SHAE_05), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme
cov_SHAE_06 <- subset(coverage, V7=="SHAE_06")
cov_SHAE_06<- cov_SHAE_06[order(cov_SHAE_06$V1, cov_SHAE_06$V2),]
cov_SHAE_06$ID <- seq.int(nrow(cov_SHAE_06))

axisdf = cov_SHAE_06 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
SHAE_06.out <- ggplot(subset(cov_SHAE_06), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("grey75", "grey40"),8)) +
  scale_x_continuous( label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) +
  theme_bw() +
  selection_theme

axisdf = cov_SHAE_06 %>% group_by(V1) %>% summarize(center=( max(ID) + min(ID) ) / 2 )
axis <- ggplot(subset(cov_SHAE_06), aes(x=ID, y=(V4))) +
  geom_point( aes(color=as.factor(V1)),size=0.00001) +
  ylab("Coverage") +
  scale_color_manual(values = rep(c("white", "white"),8)) +
  scale_x_continuous(label = axisdf$V1, breaks= axisdf$center, expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,200)) + theme_nothing()+
  theme(axis.text.x=element_text(size=6, color="black", face="bold"), 
        axis.text.y=element_text(size=5, color="black", face="bold"),
        axis.title.y=element_text(size=5, face="bold"))
```
## Merge plots
```
a22 <- plot_grid(SBOV_01.out,SBOV_02.out,
          SBOV_01.out,SBOV_02.out,
          SCUR_01.out,SGUI_01.out,
          SHAE_01.out,SHAE_02.out,
          SHAE_06.out,SHAE_04.out,
          SHAE_05.out,SINT_01.out,
          SMAT_01.out,SMAG_01.out,
          BK16_B3_01.out,BK16_B3_02.out,
          BK16_B8_01.out,BK16_B8_02.out,
          BK16_B8_03.out,BK16_B8_04.out, 
          BK16_B8_05.out, BK16_B8_06.out,axis,axis,
          ncol=2,labels=c("","",
                          "SBOV_01","SBOV_02",
                   "SCUR_01","SGUI_01",
                   "SHAE_01","SHAE_02",
                   "SHAE_03","SHAE_04",
                   "SHAE_05","SINT_01",
                   "SMAT_01","SMAG_01",
                   "BK16_B3_01","BK16_B3_02",
                   "BK16_B8_01","BK16_B8_02",
                   "BK16_B8_03","BK16_B8_04",
                   "BK16_B8_05","BK16_B8_06","",""), 
          scale = 0.8,label_size = 8, vjust=-0.125, align = "v")

a23<-plot_grid(RT15_B6_01.out,BK16_B8_07.out,
          RT15_B6_01.out,BK16_B8_05.out, 
          BK16_B8_06.out, BK16_B8_07.out,
          BK16_B8_08.out,BK16_B8_09.out,
          BK16_B8_10.out,BK16_B8_11.out,
          BK16_B8_12.out,BK16_B8_13.out,
          BK16_B8_14.out,BK16_B8_15.out, 
          BK16_B8_17.out, BK16_B8_19.out, 
          RT15_B3_01.out,RT15_B3_01.out,
          BK16_B8_14.out,BK16_B8_15.out, 
          BK16_B8_17.out, BK16_B8_19.out,
          axis,axis,
          ncol=2,labels=c("","",
                          "RT15_B6_01","BK16_B8_05",
                          "BK16_B8_06","BK16_B8_07",
                          "BK16_B8_08","BK16_B8_09",
                          "BK16_B8_10","BK16_B8_11",
                          "BK16_B8_12","BK16_B8_13",
                          "BK16_B8_14","BK16_B8_15",
                          "BK16_B8_17","BK16_B8_19",
                          "RT15_B3_01","",
                          "BK16_B8_14","BK16_B8_15",
                          "BK16_B8_17","BK16_B8_19","",""), 
          scale = 0.8,label_size = 8, vjust=-0.125, align = "v")
```
