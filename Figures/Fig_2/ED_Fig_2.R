library(tidyverse)

df_residuals <- read.table("/MyPath/df_residuals_20260604_50kb.bed",
                           sep="\t", header=TRUE)
df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]

df_residuals <- df_residuals %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),
                  ordered = TRUE))



df_residuals$seg_1 <- -1
df_residuals <- cbind(df_residuals, replicate(8, df_residuals$seg_1))
n <- ncol(df_residuals)
names(df_residuals)[(n-7):n] <- c("Speckle_seg", "Interior_Act1_seg", "Interior_Act2_seg", "Interior_Repr1_seg","Interior_Repr2_seg","Near_Lm1_seg","Near_Lm2_seg","Lamina_seg")


df_residuals$Speckle_seg <- ifelse(df_residuals$SPIN == "Speckle",df_residuals$Speckle_seg, 10)
df_residuals$Interior_Act1_seg <- ifelse(df_residuals$SPIN == "Interior_Act1",df_residuals$Interior_Act1_seg, 10)
df_residuals$Interior_Act2_seg <- ifelse(df_residuals$SPIN == "Interior_Act2",df_residuals$Interior_Act2_seg, 10)
df_residuals$Interior_Repr1_seg <- ifelse(df_residuals$SPIN == "Interior_Repr1",df_residuals$Interior_Repr1_seg, 10)
df_residuals$Interior_Repr2_seg <- ifelse(df_residuals$SPIN == "Interior_Repr2",df_residuals$Interior_Repr2_seg, 10)
df_residuals$Near_Lm1_seg <- ifelse(df_residuals$SPIN == "Near_Lm1",df_residuals$Near_Lm1_seg, 10)
df_residuals$Near_Lm2_seg <- ifelse(df_residuals$SPIN == "Near_Lm2",df_residuals$Near_Lm2_seg, 10)
df_residuals$Lamina_seg <- ifelse(df_residuals$SPIN == "Lamina",df_residuals$Lamina_seg, 10)


df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]



###############################################################################################################################################################################################################################
######## ED Fig.2b was generated using 'Fig_2.R'
######## ED Fig.2c-d were generated using 'plot_bioanalyzer_traces.R'
######## ED Fig. 2e and 2g -  Global violins    ########################

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif','LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif','LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628_dif')]


NTTD <- NT[,c(4, 6)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "TPL/DRB"
NTTD <- NTTD[NTTD$SPIN ==  c("Speckle", "Interior_Act1", "Interior_Act2","Interior_Repr1", "Interior_Repr2","Near_Lm1", "Near_Lm2", "Lamina", "Lamina_Like"),]

pdf("LOSr_allSPIN_global_NTTD1_violin_alt.pdf", width=10, height=5)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2","grey")) +
  ylim(-0.215, 0.12) +
  geom_hline(yintercept = 0, color="black", linetype=2, size=1.2, lineend = "round") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size = 20, margin = margin(r = 8))) +
  theme(axis.line   = element_line(linewidth = 2, lineend = "round"),
        axis.ticks  = element_line(linewidth = 2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

aggregate(LOS_residuals ~ SPIN, data = na.omit(NTTD), FUN = mean)



NTTD <- NT[,c(4, 7)]
colnames(NTTD)[2] <- "LOS_residuals"

#Violin plotting the LOS residuals for dTAG 6h
pdf("LOSr_allSPIN_global_U1AMO4h_violin_alt.pdf", width=10, height=5)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2","grey")) +
  ylim(-0.215, 0.12) +
  geom_hline(yintercept = 0, color="black", linetype=2, size=1.2, lineend = "round") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size = 20, margin = margin(r = 8))) +
  theme(axis.line   = element_line(linewidth = 2, lineend = "round"),
        axis.ticks  = element_line(linewidth = 2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

aggregate(LOS_residuals ~ SPIN, data = na.omit(NTTD), FUN = mean)

##########




####### ED Fig.2f
######### Crane maps were generated by "/sharehome/dlafonta/manuscripts/RNA/Figure_RNAtreat/heatmap/temp/Fig2_CraneHeatmapObsExp.ipynb"
############## Cranemaps were generated by "/sharehome/dlafonta/manuscripts/RNA/Ext_data_Fig_Hi-C/chromosome_CraneHeatmapObsExp_TreatmentsMock_LC_Hi-C.ipynb" and "/sharehome/dlafonta/manuscripts/RNA/Ext_data_Fig_Hi-C/chromosome_CraneHeatmapObsExp_with_tweaking_20200103_LC_Hi-C_NTTD_M_diff.ipynb"

seg_size <- 30
u <- 2.4
lw <- 4.2

x <-df_residuals[df_residuals$chrom == "chr11", ]
pdf("PC1_SPIN_segs_track_ch11_xaxis_alt.pdf", height = 4.2, width = 9)
par(mfrow=c(2,1), mar=c(0, 4, 0, 2) + 0.1, oma=c(5, 0, 0, 0),mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1.5,5))

Speckle <- x$Speckle_seg
segments(x0 = x$end/1000000,                   
         y0 = Speckle,
         x1 = x$start/1000000,
         y1 = Speckle,
         lend=1,lwd=seg_size, col ="#9E0142")

Interior_Act1 <- x$Interior_Act1_seg
segments(x0 = x$end/1000000,                   
         y0 = Interior_Act1,
         x1 = x$start/1000000,
         y1 = Interior_Act1,
         lend=1,lwd=seg_size, col="#F46D43")

Interior_Act2 <- x$Interior_Act2_seg
segments(x0 = x$end/1000000,                   
         y0 = Interior_Act2,
         x1 = x$start/1000000,
         y1 = Interior_Act2,
         lend=1,lwd=seg_size, col="#FDAE61")

Interior_Repr1 <- x$Interior_Repr1_seg
segments(x0 = x$end/1000000,                   
         y0 = Interior_Repr1,
         x1 = x$start/1000000,
         y1 = Interior_Repr1,
         lend=1,lwd=seg_size, col="#FFFFBF")

Interior_Repr2 <- x$Interior_Repr2_seg
segments(x0 = x$end/1000000,                   
         y0 = Interior_Repr2,
         x1 = x$start/1000000,
         y1 = Interior_Repr2,
         lend=1,lwd=seg_size, col="#ABDDA4")

Near_Lm1 <- x$Near_Lm1_seg
segments(x0 = x$end/1000000,                   
         y0 = Near_Lm1,
         x1 = x$start/1000000,
         y1 = Near_Lm1,
         lend=1,lwd=seg_size, col="#66C2A5")

Near_Lm2 <- x$Near_Lm2_seg
segments(x0 = x$end/1000000,                   
         y0 = Near_Lm2,
         x1 = x$start/1000000,
         y1 = Near_Lm2,
         lend=1,lwd=seg_size, col="#3288BD")

Lamina <- x$Lamina_seg
segments(x0 = x$end/1000000,                   
         y0 = Lamina,
         x1 = x$start/1000000,
         y1 = Lamina,
         lend=1,lwd=seg_size, col="#5E4FA2")

plot(x$start/1000000, x$PC1, type="n",
     xlab = "", ylab= "", axes=FALSE)
A <- x$PC1
A[is.na(A)] <- 0
A[A<0] <- 0
polygon(c(x$start/1000000, rev(x$start/1000000)), 
        c(A, rep(0, length(A))), col="red",
        border=NA)
B <- x$PC1
B[is.na(B)] <- 0
B[B>0] <- 0
polygon(c(x$start/1000000, rev(x$start/1000000)), 
        c(B, rep(0, length(B))), col="blue",
        border=NA)
axis(1, lwd=lw, cex.axis=2.2, labels=FALSE) 
axis(2, lwd=lw, cex.axis=u)
box(bty="l", lwd=lw)


dev.off()



x <- read.table("/MyPath/LOS_residuals_range2Mb_K562_contAMO120mPD_DpnII_R4_20240828.bedGraph",
                sep="\t", header=TRUE)
x <-x[x$chrom == "chr11", ]
pdf("LOSres_contAMO_track_ch11_xaxis_150kb.pdf", height = 2.4, width = 9)
par(mfrow=c(1,1), mar=c(2, 4, 2, 2) + 0.1)

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_contAMO120mPD_DpnII_R4_20240828, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=5, cex.axis=2.2, labels=FALSE) 
axis(2, lwd=5, cex.axis=u)
box(bty="l", lwd=5)

dev.off()


y <- read.table("/MyPath/LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828.bedGraph",
                sep="\t", header=TRUE)
y <-y[y$chrom == "chr11", ]
pdf("LOSres_U1AMO_track_ch11_xaxis_150kb.pdf", height = 2.4, width = 9)
par(mfrow=c(1,1), mar=c(2, 4, 2, 2) + 0.1)

plot(y$start/1000000, y$LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=5, cex.axis=2.2, labels=FALSE) 
axis(2, lwd=5, cex.axis=u)
box(bty="l", lwd=5)

dev.off()



x <- read.table("/MyPath/LOS_residuals_range2Mb_K562_contAMO120mPD_DpnII_R4_20240828.bedGraph",
                sep="\t", header=TRUE)
x <-x[x$chrom == "chr11", ]

y <- read.table("/MyPath/LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828.bedGraph",
                sep="\t", header=TRUE)
y <-y[y$chrom == "chr11", ]

z <- merge(x,y, by = c("chrom","start","end"),all.x= TRUE)
z$dif <- z$LOS_residuals_range2Mb_K562_contAMO120mPD_DpnII_R4_20240828 - z$LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828
z$chrom <- factor(z$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
z <- z[order(z$chrom, z$start,z$end),]
rownames(z) <- NULL

pdf("LOSres_U1AMOdif_track_ch11_xaxis_150kb.pdf", height = 2.5, width = 9)
par(mfrow=c(1,1), mar=c(2.5, 4, 2, 2) + 0.1)

plot(z$start/1000000, z$dif, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=5, cex.axis=2.2, labels=TRUE,mgp=c(3, 1.5, 0)) 
axis(2, lwd=5, cex.axis=u)
box(bty="l", lwd=5)


dev.off()

############################

