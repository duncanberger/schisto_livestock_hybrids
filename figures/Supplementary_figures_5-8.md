# Supplementary figures 5-8
## Load libraries
```
library("ggplot2")
library("reshape2")
library("dplyr")
```
## Make plots
```
# Load heterozygous site AD values
hets <- read.table("AD.all.txt")

# Load a list of samples (to loop over)
sl <- read.table("sample.hete.list")
sl2 <- as.character(unlist(sl, use.names=FALSE))

# Plot for each sample (using internal sample ID's)
for (i in seq_along(sl2)) {
    A <- ggplot(data=subset(hets, V1==sl2[i])) + 
      geom_hex(aes(y=(V4),x=(V4+V3),colour = ..count..),bins=70,) + 
      theme_bw() +
      xlab("Depth") +
      xlim(0,150) +
      ylim(0,150) +
      ylab("Depth (minor allele)") +
      geom_abline(slope = 1) +
      geom_abline(slope = 0.4) +
      scale_fill_viridis() +
      scale_color_viridis() +
      theme(panel.grid=element_blank(), 
            axis.title=element_text(face="bold", size=6, color="black"),
            legend.text = element_text(face="bold"), 
            axis.title.x = element_text(face="bold", size=6, color="black"),
            axis.title.y = element_text(face="bold", size=6, color="black"),
            legend.position = "none", 
            axis.text = element_text(face = "bold", size=5, color="black"))
    B <- ggplot(data=subset(hets, V1==sl2[i])) + 
      geom_histogram(aes(x=(V4/(V4+V3))),binwidth=0.01) + 
      theme_bw() +
      xlab("Proportion of minor allele reads") +
      ylab("Frequency") +
      theme(panel.grid=element_blank(), 
            axis.title=element_text(face="bold", size=6, color="black"),
            legend.text = element_text(face="bold"), 
            axis.title.x = element_text(face="bold", size=6, color="black"),
            axis.title.y = element_text(face="bold", size=6, color="black"),
            legend.position = "none", 
            axis.text = element_text(face = "bold", size=5, color="black")) 
  test_grid <- plot_grid(A,B, nrow=1, align= "v")
  assign(sprintf("qc_%s", sl2[i]), test_grid)
}
```
## Merge plots 
```
g1 <- plot_grid(qc_3692Sbov229274_1,qc_ERR3061640,qc_2525Scur229274_2,qc_Sgui2657_1_236593,
          qc_SRR8284792,qc_SRR8284795,qc_SRR8284796,qc_SRR8284797,qc_Sh100XfemaleNoPCR,
          qc_Sint4003_236593,qc_Smag2991_3_236593,`qc_Smat-2767PC_234514`, 
          labels = c("SBOV_01","SBOV_02","SCUR_01","SGUI_01",
                     "SHAE_01","SHAE_02","SHAE_03","SHAE_04","SHAE_05",
                     "SINT_01","SMAG_01","SMAT_01"),label_size = 8,scale=0.85, label_y = 1.0, ncol=2)

g2 <- plot_grid(qc_hybrid_schisto7118176,qc_hybrid_schisto7118177,qc_hybrid_schisto7118178,qc_hybrid_schisto7118179,
          qc_hybrid_schisto7118180,qc_hybrid_schisto7118181,qc_hybrid_schisto7118182,qc_hybrid_schisto7118183,
          qc_hybrid_schisto7714151,qc_hybrid_schisto7714152,qc_hybrid_schisto7714153,qc_hybrid_schisto7714154,
          ncol=2,
          labels = c("BK16_B3_01","BK16_B3_02","BK16_B8_01","BK16_B8_02",
                     "BK16_B8_03","BK16_B8_04","BK16_B8_05","BK16_B8_06",
                     "BK16_B8_07","BK16_B8_08","BK16_B8_09","BK16_B8_10"),
          label_size = 8,scale=0.85, label_y = 1.0)

g3 <- plot_grid(qc_hybrid_schisto7714155,qc_hybrid_schisto7714156,qc_hybrid_schisto7714157,qc_hybrid_schisto7714158,
          qc_hybrid_schisto7714159,qc_hybrid_schisto7714160,qc_hybrid_schisto7714161,qc_hybrid_schisto7714163,
          qc_hybrid_schisto7717246,qc_hybrid_schisto7717247,qc_hybrid_schisto7148988,qc_hybrid_schisto7148989,
          ncol=2, labels = c("BK16_B8_11","BK16_B8_12","BK16_B8_13","BK16_B8_14","BK16_B8_15",
                             "BK16_B8_16 (excluded)","BK16_B8_17","BK16_B8_19","RT15_B3_01",
                             "RT15_B3_02 (excluded)","RT15_B6_01","RT15_B6_02 (excluded)"),
          label_size = 8,scale=0.85, label_y = 1.0)
