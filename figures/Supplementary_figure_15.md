# Supplementary figure 15
## Load libraries
```
library("junctions")
library("ggplot2")
```
## For each level of initial heterozygosity (H_0), and number of junctions (j) 
```
for(i in 1:100) {
  x <- estimate_time(J = i, N = Inf, R = Inf, H_0 = 0.75, C = 9.8)
  print(x)
}

# Write the results to file - jcounts.csv
```
## Plot
```
jcurve <- read.csv("jcounts.csv", header=TRUE, sep='\t')

junctions <- ggplot(data=jcurve) + 
  geom_line(aes(x=junctions, y=generations, color=as.factor(Het), fill=as.factor(Het)), size=0.5) +
  scale_x_continuous(limits=c(0,30), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,20)) + 
  xlab("Number of junctions") + ylab("Number of generations") + labs(colour = "Initial heterozygosity") +
  theme_bw() + theme(axis.title=element_text(face="bold",size=7, color="black"),
                     axis.text=element_text(face="bold",size=6, color="black"),
                     panel.background=element_blank(),
                     panel.border = element_rect(color="#4c4c4c",fill=NA),
                     panel.grid=element_blank(),
                     legend.text=element_text(size=6, face="bold"),
                     legend.title=element_text(face="bold", size=7,color="black"))
```
