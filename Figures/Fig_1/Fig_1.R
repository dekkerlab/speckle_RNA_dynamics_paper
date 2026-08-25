library(tidyverse)
library(ggplot2)
library(Hmisc)




df_residuals <- read.table("/MyPath/",
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
df_residuals <- cbind(df_residuals, replicate(9, df_residuals$seg_1))
n <- ncol(df_residuals)
names(df_residuals)[(n-8):n] <- c("Speckle_seg", "Interior_Act1_seg", "Interior_Act2_seg", "Interior_Repr1_seg","Interior_Repr2_seg","Near_Lm1_seg","Near_Lm2_seg","Lamina_seg","Lamina_Like_seg")
 
 
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


get_moving_average <- function(x, y, n, s, dec) {
  ## Args:
  #  x: unordered numeric vector to slide windows across
  #  y: unordered numeric values to mean over x window
  #  n: window size
  #  s: step size
  #  dec: number of decimal places to round for merge
  # Returns:
  # df: dataframe with columns:
  #   w: window centers
  #   mu: moving average
  if ((sum(is.na(x)) > 0) | (sum(is.na(y)) > 0)) {
    stop("Remove NAs")
  }
  if (length(x)!=length(y)) {
    stop("Unequal vector lengths")
  }
  xmin <- round(min(x, na.rm=TRUE), dec)
  xmax <- round(max(x, na.rm=TRUE), dec)
  w <- seq(xmin, xmax, by=s)
  mu <- c()
  for (i in w) {
    m <- mean(y[x > (i-(n/2)) & x < (i + (n/2))])
    mu <- c(mu, m)
  }
  df <- data.frame(w, mu)
  return(df)
}

get_ma_residuals <- function(o, df, dec) {
  # Get moving average residuals
  ## Args:
  #  o: original dataframe with cols:
  #   column1: chrom
  #   column2: start
  #   column3: end
  #   column4: numeric vector original x-axis
  #   column5:  numeric vector original y-axis
  #  df: dataframe output by get_moving_average
  #      function
  #  dec: number of decimal places to round for merge
  # Returns:
  #  r: dataframe of residuals and locations
  colnames(o) <- c("chrom", "start", "end", "w", "y")
  o["w"] <- round(o["w"], dec)
  df["w"] <- round(df["w"], dec)
  m <- merge(o, df, by="w")
  m_resid <- m$y-m$mu
  r <- cbind(m[c("chrom", "start", "end")], m_resid)
  return(r)
}


#####Figure 1#####
## Fig. 1 - all heatmaps were generated using "/sharehome/dlafonta/manuscripts/RNA/Figure_sp_stability/heatmap/Fig1_heatmaps.ipynb"
## Fig.1a  - - PC1-colored track (chr11) 
seg_size <- 20


setwd("/MyPath/")
x <-df_residuals[df_residuals$chrom == "chr11", ]
pdf("PC1_SPIN_segs_track_ch11_xaxis_alt.pdf", height = 3.6, width = 9)
par(mfrow=c(2,1), mar=c(0, 4, 0, 2) + 0.1, oma=c(5, 0, 0, 0),mgp=c(3, 1.5, 0))
plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1.5,4))

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
axis(1, lwd=5, cex.axis=2, labels=FALSE) 
axis(2, lwd=5, cex.axis=2.2)

dev.off()


## Fig.1a - Matrix x-axis
setwd("/MyPath/")
x <-df_residuals[df_residuals$chrom == "chr11", ]
pdf(paste("xAxis_track for matrix.pdf", sep = ""), height = 1.9, width = 9)
par(mfrow=c(1,1), mar=c(3, 4, 1, 2) + 0.1, mgp=c(3, 1.5, 0))
plot(x$start/1000000, x$LOS_range2Mb_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(0, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=5, cex.axis=2.5, labels=TRUE) 

dev.off()  
  
  
## Fig.1b - PC1-colored track (chr11)
x <-df_residuals[df_residuals$chrom == "chr11", ]
pdf("PC1_SPIN_segs_track_ch11_xaxis1b_alt.pdf", height = 3.6, width = 9)
par(mfrow=c(2,1), mar=c(0, 4, 0, 2) + 0.1, oma=c(5, 0, 0, 0),mgp=c(3, 1.5, 0))
plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1.5,4))

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
axis(1, lwd=5, cex.axis=2, labels=FALSE) 

dev.off()
  
#####
  

  
## Fig.1b - Raw LOS (chr11) 
setwd("/MyPath/")
x <-df_residuals[df_residuals$chrom == "chr11", ]
pdf(paste("Raw_LOS_120m_20250326_chr11_track.pdf", sep = ""), height = 2.2, width = 9)
par(mfrow=c(1,1), mar=c(3, 4, 1, 2) + 0.1, mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_range2Mb_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(0.58, 0.88), xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=5, cex.axis=2.5, labels=TRUE) 
axis(2, lwd=5, cex.axis=2.2)

dev.off()

## Fig.1c - Raw LOS HS (chr11 zoom in)
setwd("/MyPath/")
x <-df_residuals[df_residuals$chrom == "chr11", ]
x <- x[x$start>55000000 & x$end<85000000,]
pdf(paste("Raw_LOS_120m_20250326_chr11_track_zoom.pdf", sep = ""), height = 2.2, width = 9)
par(mfrow=c(1,1), mar=c(3, 4, 1, 2) + 0.1, mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_range2Mb_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(0.58, 0.88), xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=5, cex.axis=2.5, labels=TRUE) 

dev.off()


### Fig.1c - Segments and PC1 (zoom)

seg_size <- 21.5

x <-df_residuals[df_residuals$chrom == "chr11", ]
x <- x[x$start>55000000 & x$end<85000000,]
pdf("IPG_segs_PC1_track_ch11_xaxis_zoom_alt.pdf", height = 3.65, width = 9.4)
par(mfrow=c(2,1), mar=c(0.5, 4, 0, 2) + 0.1, oma=c(5, 0, 0, 0),mgp=c(3, 1.5, 0))
plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1.5,4))

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
     ylim=c(-1,1.5), xlab = "", ylab= "", axes=FALSE)
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
axis(1, lwd=5, cex.axis=2, labels=FALSE) 
axis(2, lwd=5, cex.axis=2, labels=TRUE) 
dev.off()


#######################


#################   Scatter plots #####################################
setwd("/MyPath/")
NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','signal_K562_2hPD_DpnII_R2_20250326','LOS_range2Mb_K562_2hPD_DpnII_R2_20250326')]
NT$PC_sign <- ifelse(NT$PC1 > 0, "A","B")
NT <- NT[!is.na(NT[,colnames(NT)[6]]),]
NT <- NT[!is.na(NT$SPIN),]
NT <- NT[!is.na(NT[,colnames(NT)[7]]),]


### Fig.1d - Scatter plot of LOS vs DpnII-seq signal colored by SPIN 

n <- mean(na.omit(NT[,colnames(NT)[6]]))/4
ma <- get_moving_average(NT[,colnames(NT)[6]], NT[,colnames(NT)[7]], n, 1, 0)
los_r <- get_ma_residuals(NT[c("chrom", "start", "end", colnames(NT)[6],colnames(NT)[6])], ma, 0)
pdf(paste("Signal_vs_",colnames(NT)[7],"_altLL.pdf", sep = ""), width=3, height=3)
ggplot() + 
  geom_point(data = NT, aes(x=NT[,colnames(NT)[6]], y=NT[,colnames(NT)[7]], color=SPIN),size=0.5) +
  scale_color_manual(values = c("Speckle" = "#9E0142", "Interior_Act1" ="#F46D43", "Interior_Act2" ="#FDAE61", #"Interior_Act3"="#FDAE61", 
                                "Interior_Repr1"="#FFFFBF", "Interior_Repr2"="#ABDDA4", "Near_Lm1"="#66C2A5", "Near_Lm2"="#3288BD", "Lamina"="#5E4FA2")) +
  labs(x="", y = "") +
  geom_line(data= ma, aes(x= w, y=mu)) +
  lims(x= c(0,2700), y = c(0.45, 0.95))+
  theme_classic() +
  theme(axis.text=element_text(size=13))+
  theme(axis.line = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks.length = unit(0.2, "cm")) +
  theme(legend.position = 'none')
dev.off()


##Fig.1e  - Violin plot of LOS residuals subset by PC1 sign 

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326')]
colnames(NT)[6] <- "LOS_residuals"
NT$PC_sign <- ifelse(NT$PC1 > 0, "A","B")

setwd("/MyPath/")
pdf("LOSres_2h_20250326_by_PC1_violin.pdf", width=2.4, height=2.5)
ggplot(na.omit(NT), aes(x=PC_sign, y=LOS_residuals, fill=PC_sign)) + 
  geom_violin(trim=FALSE, show.legend = FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values = c("red","blue")) +
  ylim(-0.33, 0.27) +
  labs(x="Compartment", y = "LOS")+
  geom_hline(yintercept = 0,color="black",linetype=2) +
  guides(fill=guide_legend(title="Compartment"))+
  theme_classic()+
  theme(axis.text=element_text(size=18))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.line = element_line(linewidth = 1, lineend = "round"), 
        axis.ticks = element_line(linewidth = 1, lineend = "round"),
        axis.ticks.length = unit(0.2, "cm"))

dev.off()


### Fig.1f - Violin plotting the LOS residuals by SPIN state
setwd("/MyPath/")
pdf("LOSres_2h_20250326_by_SPIN_altLL.pdf", width=5.8, height=2.5)
ggplot(na.omit(NT), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.33, 0.27) +
  labs(x="SPIN state", y = "Loss Of Structure (LOS)")+
  geom_hline(yintercept = 0,color="black",linetype=2) +
  theme_classic()+
  theme(axis.text=element_text(size=18))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.line = element_line(linewidth = 0.7, lineend = "round"),
        axis.ticks = element_line(linewidth = 0.7, lineend = "round"),
        axis.ticks.length = unit(0.2, "cm"))

dev.off()


###### Calculating % genome and % in A

spin_order <- c("Speckle", "Interior_Act1", "Interior_Act2",
                "Interior_Repr1", "Interior_Repr2", "Near_Lm1", "Near_Lm2", "Lamina")

spin_counts <- df_residuals %>%
  filter(!is.na(SPIN)) %>%
  group_by(SPIN) %>%
  summarise(
    total = n(),
    PC1_positive = sum(PC1 > 0, na.rm = TRUE),
  ) %>%
  mutate(
    pct_positive = round(PC1_positive / total * 100, 1),
    SPIN = factor(SPIN, levels = spin_order)
  ) %>%
  arrange(SPIN)

spin_counts$genome_pct <- spin_counts$total/sum(spin_counts$total)


spin_counts <- df_residuals %>%
  filter(!is.na(SPIN)) %>%
  group_by(SPIN) %>%
  summarise(
    total = n(),
    PC1_positive = sum(PC1 > 0, na.rm = TRUE),
    PC1_negative = sum(PC1 < 0, na.rm = TRUE),
    PC1_zero_or_NA = sum(PC1 == 0 | is.na(PC1))
  ) %>%
  mutate(
    pct_positive = round(PC1_positive / total * 100, 1),
    pct_negative = round(PC1_negative / total * 100, 1),
    SPIN = factor(SPIN, levels = spin_order)
  ) %>%
  arrange(SPIN)


### Fig.1g - SPIN, PC1, LOS, DpnII-seq and LOS residuals tracks
setwd("/MyPath/")
x <-df_residuals[df_residuals$chrom == "chr11", ]
pdf("LOS_correction_range2Mb_tracks_20250326_2h_chr11_alt.pdf", height = 5.5, width = 8)
par(mfrow=c(5,1), mar=c(2, 4, 0.5, 2) + 0.1)

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
     xlab = "", ylab= "", axes=FALSE, cex.lab=1.2)
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
axis(1, lwd=2, cex.axis=1.8, labels=FALSE) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l", lwd=2)

plot(x$start/1000000, x$LOS_range2Mb_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(0.58, 0.87), xlab = "", ylab= "", axes=FALSE, lwd=0.5, cex.lab=1.2)
axis(1, lwd=2, cex.axis=1.8, labels=FALSE) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l", lwd=2)

plot(x$start/1000000, x$signal_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(0, 2000), xlab = "", ylab= "", axes=FALSE, lwd=0.5, cex.lab=1.2)
axis(1, lwd=2, cex.axis=1.8, labels=FALSE) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l", lwd=2)

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(-0.22, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=0.5, cex.lab=1.2)
axis(1, lwd=2, cex.axis=1.8, labels=TRUE) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l", lwd=2)


dev.off()


###

###########
