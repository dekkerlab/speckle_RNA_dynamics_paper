library(tidyverse)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggbeeswarm)
library(rtracklayer)
library(data.table)
library(readxl)

df_residuals <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/dataframes/df_residuals_20260604_50kb.bed",
                           sep="\t", header=TRUE)
df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]


#df_residuals$SPIN <- factor(df_residuals$SPIN, order = TRUE, levels =c("Speckle", "Interior_Act1", "Interior_Act2", "Interior_Act3", "Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina","Lamina_like"))
df_residuals <- df_residuals %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina","Lamina_like"),
                  ordered = TRUE))

df_residuals$seg_1 <- -1
df_residuals <- cbind(df_residuals, replicate(9, df_residuals$seg_1))
n <- ncol(df_residuals)
names(df_residuals)[(n-8):n] <- c("Speckle_seg", "Interior_Act1_seg", "Interior_Act2_seg", "Interior_Act3_seg", "Interior_Repr1_seg","Interior_Repr2_seg","Near_Lm1_seg","Near_Lm2_seg","Lamina_seg")


df_residuals$Speckle_seg <- ifelse(df_residuals$SPIN == "Speckle",df_residuals$Speckle_seg, 10)
df_residuals$Interior_Act1_seg <- ifelse(df_residuals$SPIN == "Interior_Act1",df_residuals$Interior_Act1_seg, 10)
df_residuals$Interior_Act2_seg <- ifelse(df_residuals$SPIN == "Interior_Act2",df_residuals$Interior_Act2_seg, 10)
df_residuals$Interior_Act3_seg <- ifelse(df_residuals$SPIN == "Interior_Act3",df_residuals$Interior_Act3_seg, 10)
df_residuals$Interior_Repr1_seg <- ifelse(df_residuals$SPIN == "Interior_Repr1",df_residuals$Interior_Repr1_seg, 10)
df_residuals$Interior_Repr2_seg <- ifelse(df_residuals$SPIN == "Interior_Repr2",df_residuals$Interior_Repr2_seg, 10)
df_residuals$Near_Lm1_seg <- ifelse(df_residuals$SPIN == "Near_Lm1",df_residuals$Near_Lm1_seg, 10)
df_residuals$Near_Lm2_seg <- ifelse(df_residuals$SPIN == "Near_Lm2",df_residuals$Near_Lm2_seg, 10)
df_residuals$Lamina_seg <- ifelse(df_residuals$SPIN == "Lamina",df_residuals$Lamina_seg, 10)


CPM <- import.bw("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/SLAMseq/coverage/HS37_20250528_50000.bw")

# Convert to data frame (bed-like)
CPM <- as.data.frame(CPM)[, c("seqnames","start","end","score")]
colnames(CPM)[1] <- "chrom"
colnames(CPM)[4] <- "CPM"
CPM$chrom <- paste0("chr", CPM$chrom)
CPM$start <- CPM$start -1


df_residuals <- merge(df_residuals,CPM, by= c("chrom","start","end"), all.x = TRUE)

###### Count genes per bin

genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
bins <- fread("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Sequencing/data/binning/hg38_50kb.bed", col.names = c("chrom","start","end"))
bins_gr <- makeGRangesFromDataFrame(bins, starts.in.df.are.0based = TRUE)
bins$gene_count <- countOverlaps(bins_gr, genes_gr)



df_residuals <- merge(df_residuals,bins, by= c("chrom","start","end"), all.x = TRUE)

df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]




########################  Tracks for regions targeted by DNA-FISH ######################

##### Extended Data Fig. 3a

#### Tracks for chr19:45M

i <- "chr3"
start <- 45000000
end <- 60000000

#### sites for chr19:45M
site1 <- 49682567
site2 <- 49722567

seg_size <- 30

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/plots/")

x <-df_residuals[df_residuals$chrom == i, ]
x <- x[x$start>=start & x$end<=end,]

lw <- 3.7 
ft <- 2.6
  
pdf(paste("DNA-FISH_track_",i,"_",site1,"_alt.pdf",sep=""), height = 8, width = 7.5)
par(mfrow=c(5,1), mar=c(2.5, 4, 0.5, 2) + 0.1,mgp=c(3, 1.5, 0))
plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1,2))

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
     ylim=c(-1.5, 1.5), xlab = "", ylab= "", axes=FALSE)
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
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)


plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif, type="l", col="black", 
     ylim=c(-0.15, 0.05), xlab = "", ylab= "", axes=FALSE, lwd=1.5, yaxt='n')
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft, at=c(-0.15, -0.10, -0.05, 0, 0.05), labels=c("-0.15", "", "", "0", ""))
abline(h = 0,col = "black", lty=2)
abline(v = site1/1000000,col = "red", lty="dashed")
abline(v = site2/1000000,col = "red", lty="dashed")
box(bty="l", lwd=lw)

plot(x$start/1000000, x$CPM, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

plot(x$start/1000000, x$gene_count, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=lw, cex.axis=ft, labels=TRUE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

dev.off()


########################  Tracks for regions targeted by DNA-FISH ######################
##### Extended Data Fig. 3b

#### Tracks for chr19:45M

i <- "chr15"
start <- 68000000
end <- 83000000

#### sites for chr19:45M
site1 <- 75467659
site2 <- 75507659

lw <- 3.7 
ft <- 2.6

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/plots/")

x <-df_residuals[df_residuals$chrom == i, ]
x <- x[x$start>=start & x$end<=end,]

pdf(paste("DNA-FISH_track_",i,"_",site1,"_alt.pdf",sep=""), height = 8, width = 7.5)
par(mfrow=c(5,1), mar=c(2.5, 4, 0.5, 2) + 0.1,mgp=c(3, 1.5, 0))
plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1,2))

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
     ylim=c(-1.5, 1.5), xlab = "", ylab= "", axes=FALSE)
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
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)


plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif, type="l", col="black", 
     ylim=c(-0.15, 0.05), xlab = "", ylab= "", axes=FALSE, lwd=1.5, yaxt='n')
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft, at=c(-0.15, -0.10, -0.05, 0, 0.05), labels=c("-0.15", "", "", "0", ""))
abline(h = 0,col = "black", lty=2)
abline(v = site1/1000000,col = "red", lty="dashed")
abline(v = site2/1000000,col = "red", lty="dashed")
box(bty="l", lwd=lw)

plot(x$start/1000000, x$CPM, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

plot(x$start/1000000, x$gene_count, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=lw, cex.axis=ft, labels=TRUE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

dev.off()


###########################################################################################################################################
################################## Plotting all groups and channels in  a single panel #######################################################################################


###### Fig. 3b

###### Segmentation with dilated mask example in Extended data Fig 3c was generated using "plot_segmentationMask.ipynb" (Globally normalized all channels except SON) 
######

##### Images were normalized, segmented, dilated mask applied and SON signal quantified using "FISH_quant.ipynb" (Globally normalized all channels except SON)

NT_568_rawSON <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/DMSO_568/DMSO_568_foci_SON_intensity.csv",
                            sep=",", header=TRUE)

TD_568_rawSON <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/TD_568/TD_568_foci_SON_intensity.csv",
                            sep=",", header=TRUE)

NT_647_rawSON <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/DMSO_647/DMSO_647_foci_SON_intensity.csv",
                            sep=",", header=TRUE)

TD_647_rawSON <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/TD_647/TD_647_foci_SON_intensity.csv",
                            sep=",", header=TRUE)

### Make a treatment name column
NT_568_rawSON$treatment <- "DMSO"
TD_568_rawSON$treatment <- "TD"
NT_647_rawSON$treatment <- "DMSO"
TD_647_rawSON$treatment <- "TD"

### Make channel column
NT_568_rawSON$channel <- "568"
TD_568_rawSON$channel <- "568"
NT_647_rawSON$channel <- "647"
TD_647_rawSON$channel <- "647"

## Make group column 
NT_568_rawSON$group <- "DMSO_568"
TD_568_rawSON$group <- "TD_568"
NT_647_rawSON$group <- "DMSO_647"
TD_647_rawSON$group <- "TD_647"

NTTD_allChannels <- rbind(NT_568_rawSON,TD_568_rawSON,NT_647_rawSON,TD_647_rawSON)
NTTD_allChannels$group <- factor(NTTD_allChannels$group, 
                                 levels = c("DMSO_568", "TD_568", "DMSO_647", "TD_647"))


summary_df <- NTTD_allChannels %>%
  group_by(group) %>%
  summarise(
    n_cells  = n_distinct(sub("/.*$", "", collection_name)),
    n_slices = n_distinct(slice_identifier),
    n_foci   = n(),
    .groups  = "drop"
  )

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/plots")
pdf("DNAFISH_IF_NTTD_bothChannels_rawSON_dodged_KS.pdf", width = 14, height = 5)
ggplot(NTTD_allChannels, aes(x = group, y = NTTD_allChannels[,6], fill = group)) +
  geom_beeswarm(cex = 1.2, size = 2, aes(color = group)) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  ylim(0,140) +
  scale_color_manual(values = c("red", "blue", "red", "blue")) + 
  labs(x = "", y = "Mean normalized") +
  theme_classic() +
  theme(axis.ticks = element_line(colour = "black", size = 2.1, lineend = "round"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(colour = "black", size = 2.1, lineend = "round"),
        axis.text = element_text(size = 30),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

dev.off()

groups <- levels(NTTD_allChannels$group)
col_of_interest <- 6

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    x <- NTTD_allChannels[NTTD_allChannels$group == groups[i], col_of_interest]
    y <- NTTD_allChannels[NTTD_allChannels$group == groups[j], col_of_interest]
    ks <- ks.test(x, y)
    cat(sprintf("%s vs %s: D = %s, p = %s\n", groups[i], groups[j],
                format(ks$statistic, scientific = FALSE),
                format(ks$p.value, scientific = FALSE)))
  }
}




##################################################################################################################################
############################################# plot foci area ####################################################################


##### Extended Data Fig. 3d

NT_568 <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/DMSO_568/DMSO_568_combined_measurements.csv",
                     sep=",", header=TRUE)
NT_568$condition <- "NT"
NT_568$channel <- "568"

TD_568 <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/TD_568/TD_568_combined_measurements.csv",
                     sep=",", header=TRUE)
TD_568$condition <- "TD"
TD_568$channel <- "568"


NT_647 <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/DMSO_647/DMSO_647_combined_measurements.csv",
                     sep=",", header=TRUE)
NT_647$condition <- "NT"
NT_647$channel <- "647"

TD_647 <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/TD_647/TD_647_combined_measurements.csv",
                     sep=",", header=TRUE)
TD_647$condition <- "TD"
TD_647$channel <- "647"

NTTD_568 <- rbind(NT_568, TD_568)
NTTD_647 <- rbind(NT_647, TD_647)

### Make a treatment name column
NTTD_568$treatment <- NTTD_568$condition
NTTD_568$cell_id <- paste0(NTTD_568$condition, "_", sub("/.*", "",  NTTD_568$collection_name))

NTTD_647$treatment <- NTTD_647$condition
NTTD_647$cell_id <- paste0(NTTD_647$condition, "_", sub("/.*", "",  NTTD_647$collection_name))


### Catenate 568 and 647 channels
FISH_all <- rbind(NTTD_568, NTTD_647)

####### Plot foci size for both 568 and 647 channels  #####################

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/plots")
pdf("DNAFISH_IF_NTTD_ch3FocusArea2.pdf", width=6, height=5)

NT <- FISH_all[,c('treatment', 'foci_area_px', 'channel','cell_id')]
NT <- na.omit(NT)

annotation_df <- NT %>%
  group_by(channel) %>%
  summarise(n_cells = n_distinct(cell_id), n_foci  = n(), 
            y_pos   = max(NT[,2]) * 1.05) %>%
  mutate(label = paste0("#cells: ", n_cells, "\n","#foci: ", n_foci))


ggplot(NT, aes(x = channel, y = NT[,2], fill = channel)) +
  geom_beeswarm(cex = 2, size = 2, aes(color = channel)) +
  geom_boxplot(width = 0.05,fill = "white", outlier.shape = NA) +
  geom_text(data = annotation_df, aes(x = channel, y = y_pos, label = label),inherit.aes = FALSE, size = 3.5, vjust = 0) +
  scale_color_manual(values = c("red", "blue")) + 
  labs(x = "Channel", y = "focus size") +
  ggtitle(colnames(NT)[2]) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

dev.off()


NT$treatment <- as.factor(NT$treatment)
groups <- levels(NT$treatment)
col_of_interest <- 2

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    x <- NT[NT$treatment == groups[i], col_of_interest]
    y <- NT[NT$treatment == groups[j], col_of_interest]
    ks <- ks.test(x, y)
    cat(sprintf("%s vs %s: D = %s, p = %s\n", groups[i], groups[j],
                format(ks$statistic, scientific = FALSE),
                format(ks$p.value, scientific = FALSE)))
  }
}



##################################################################################################################################
############################################# plot ch2 intensity by treatment size ####################################################################


##### Extended Data Fig. 3e

NT_568 <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/DMSO_568/DMSO_568_combined_measurements.csv",
                       sep=",", header=TRUE)
NT_568$condition <- "NT"

TD_568 <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/TD_568/TD_568_combined_measurements.csv",
                sep=",", header=TRUE)
TD_568$condition <- "TD"

NTTD_568 <- rbind(NT_568,TD_568)


NTTD_568$cell_id <- sub("/.*", "",  NTTD_568$collection_name)


####### Plot foci size for both 568 and 647 channels  #####################
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/plots")
pdf("DNAFISH_IF_SONintensityByTreatment2.pdf", width=6, height=5)

NT <- NTTD_568[,c('condition', 'mean_intensity', 'cell_id')]
NT <- na.omit(NT)


annotation_df <- NT %>%
  group_by(condition) %>%
  summarise(n_cells = n_distinct(cell_id), n_foci  = n(),                     # each row = one focus
            y_pos   = max(NT[,2]) * 1.05) %>%
  mutate(label = paste0("#cells: ", n_cells, "\n","#foci: ", n_foci))


ggplot(NT, aes(x = condition, y = NT[,2], fill = condition)) +
  geom_beeswarm(cex = 2, size = 2, aes(color = condition)) +
  geom_boxplot(width = 0.05,fill = "white", outlier.shape = NA) +
  geom_text(data = annotation_df, aes(x = condition, y = y_pos, label = label),inherit.aes = FALSE, size = 3.5, vjust = 0) +
  scale_color_manual(values = c("red", "blue")) + 
  labs(x = "Channel", y = "focus size") +
  ggtitle(colnames(NT)[2]) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

dev.off()


NT$condition <- as.factor(NT$condition)
groups <- levels(NT$condition)
col_of_interest <- 2

for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    x <- NT[NT$condition == groups[i], col_of_interest]
    y <- NT[NT$condition == groups[j], col_of_interest]
    ks <- ks.test(x, y)
    cat(sprintf("%s vs %s: D = %s, p = %s\n", groups[i], groups[j],
                format(ks$statistic, scientific = FALSE),
                format(ks$p.value, scientific = FALSE)))
  }
}




##### Extended Data Fig. 3f

##################### Is TSA-seq signal at FISH loci
SONtsaseq_raw <- import.bw("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Sequencing/data/K562_data/DL/SONCond2_TSAseqV1_R1_4DNFINVM6LFT.bw")

m <- df_residuals[, c("chrom", "start", "end", "LOS_residuals_range2Mb_K562_mega120mPD_DpnII_R12_20250129")]

# keep.extra.columns = TRUE so LOS carries into gr_m
gr_bed <- makeGRangesFromDataFrame(SONtsaseq_raw , keep.extra.columns = TRUE)
gr_m   <- makeGRangesFromDataFrame(m, keep.extra.columns = TRUE)

# Find all overlapping pairs
hits    <- findOverlaps(gr_m, gr_bed)
m_idx   <- queryHits(hits)
bed_idx <- subjectHits(hits)

# Build overlaps table with coordinates, overlap length, and both value columns
overlaps <- tibble(
  m_idx   = m_idx,
  bed_idx = bed_idx
) |>
  mutate(
    chrom     = as.character(seqnames(gr_m))[m_idx],
    m_start   = start(gr_m)[m_idx],
    m_end     = end(gr_m)[m_idx],
    bed_start = start(gr_bed)[bed_idx],
    bed_end   = end(gr_bed)[bed_idx],
    .ovlp     = pmin(m_end, bed_end) - pmax(m_start, bed_start),
    score     = mcols(gr_bed)$score[bed_idx],
    LOS       = mcols(gr_m)$LOS_residuals_range2Mb_K562_mega120mPD_DpnII_R12_20250129[m_idx]
  )




result <- overlaps |>
  group_by(chrom, start = m_start, end = m_end) |>
  summarise(
    score = sum(score * .ovlp, na.rm = TRUE) / sum(.ovlp[!is.na(score)]),
    .groups = "drop"
  )

# Left join back to m to preserve any gr_m bins with no overlapping bed_df rows
result <- m |>
  left_join(result, by = c("chrom", "start", "end"))


result <- na.omit(result)




# Define regions to label
chrom1 <- "chr3"
start1 <- 49682567
end1   <- 49722567

chrom2 <- "chr15"
start2 <- 75467659
end2   <- 75507659

# Pull the matching rows from result
label_points <- result |>
  filter(
    (chrom == chrom1 & start < end1 & end > start1) |
      (chrom == chrom2 & start < end2 & end > start2)
  ) |>
  mutate(label = case_when(
    chrom == chrom1 ~ "chr3:50M",
    chrom == chrom2 ~ "chr15:75M"
  ))


setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/plots")
#Violin plotting SON TSA-seq distribution while highlithin regions that are probed by FISH
pdf("TSA-seqvsLOS_atFISHloci.pdf", width=5, height=3)
ggplot(result, aes(x = score, y = LOS_residuals_range2Mb_K562_mega120mPD_DpnII_R12_20250129)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c() +
  geom_point(data = label_points, aes(x = score, y = LOS_residuals_range2Mb_K562_mega120mPD_DpnII_R12_20250129, color = label),
             size = 2) +
  scale_color_manual(values = c("chr3:50M" = "red", "chr15:75M" = "dodgerblue")) +
  labs(x = "Score", y = "LOS residuals", color = "Region") +
  theme_classic()
dev.off()


##### Fig. 3d


########################  Tracks for regions targeted by CRISPR-SIRIUS ######################

#### Tracks for chr19:45M
seg_size <- 30
i <- "chr19"
start <- 44500000
end <- 60000000

#### sites for chr19:45M
site1 <- 45303993
site2 <- 45305361
site3 <- 49465064
site4 <- 49475276

lw <- 3
ft <- 2

setwd("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_CRISPR_SIRIUS")

x <-df_residuals[df_residuals$chrom == i, ]
x <- x[x$start>=start & x$end<=end,]

pdf(paste("CRISPR_SIRIUS_targets_tracks_",i,".pdf",sep=""), height = 6, width = 5.5)
par(mfrow=c(5,1), mar=c(3, 4, 0, 2) + 0.1)
plot(x$start/1000000, x$seg_1, type="n",
     xlab = "", ylab= "", axes=FALSE, ylim=c(-1,2))

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
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif, type="l", col="black", 
     ylim=c(-0.15, 0.05), xlab = "", ylab= "", axes=FALSE, lwd=1)
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft,  
     at=c(-0.15, -0.1, -0.05, 0, 0.05),
     labels=c(-0.15, "", "", 0, ""))
abline(h = 0,col = "black", lty=2)
abline(v = site1/1000000,col = "blue", lty=1)
abline(v = site2/1000000,col = "blue", lty=1)
abline(v = site3/1000000,col = "red", lty=1)
abline(v = site4/1000000,col = "red", lty=1)
box(bty="l", lwd=lw)

plot(x$start/1000000, x$CPM, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=2)
axis(1, lwd=lw, cex.axis=ft, labels=FALSE) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

plot(x$start/1000000, x$gene_count, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=1)
axis(1, lwd=lw, cex.axis=ft, labels=TRUE, mgp=c(3, 1.5, 0),) 
axis(2, lwd=lw, cex.axis=ft)
box(bty="l", lwd=lw)

dev.off()



#### Fig. 3f

#### Load and munge Imaging data (CRISPR-Sirius)  ###############

Rg_table <- read_excel("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/Rg-Jan10_2026_85.xlsx", skip = 3, col_names = FALSE)

# Extract the two header rows
header1 <- as.character(Rg_table[1, ])
header2 <- as.character(Rg_table[2, ])

# Combine headers into one column name
col_names <- paste(header1, header2, sep = " | ")

# Remove the two header rows and set column names
Rg_table <- Rg_table[-(1:2), ]
colnames(Rg_table) <- col_names

# Convert to long format
df_long <- Rg_table |>
  mutate(row_id = row_number()) |>
  pivot_longer(
    cols = -row_id,
    names_to = "group",
    values_to = "value"
  ) |>
  separate(group, into = c("region", "condition"), sep = " \\| ") |>
  mutate(value = as.numeric(value))

df_long <- df_long |>
  mutate(condition = case_when(
    str_detect(condition, "DMSO") ~ "DMSO",
    str_detect(condition, "TX")   ~ "TX"
  ))




Deff_table <- read_excel("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/data/LDiffB2 const_Jan10_2026.xlsx", skip = 3, col_names = FALSE)

# Extract the two header rows
header1 <- as.character(Deff_table[1, ])
header2 <- as.character(Deff_table[2, ])

# Combine headers into one column name
col_names <- paste(header1, header2, sep = " | ")

# Remove the two header rows and set column names
Deff_table <- Deff_table[-(1:2), ]
colnames(Deff_table) <- col_names

# Convert to long format
Deff_table_long <- Deff_table |>
  mutate(row_id = row_number()) |>
  pivot_longer(
    cols = -row_id,
    names_to = "group",
    values_to = "value"
  ) |>
  separate(group, into = c("region", "condition"), sep = " \\| ") |>
  mutate(value = as.numeric(value))

Deff_table_long <- Deff_table_long |>
  mutate(condition = case_when(
    str_detect(condition, "DMSO") ~ "DMSO",
    str_detect(condition, "TX")   ~ "TX"
  ))




#### Figure 4e - Radius of gyration boxplots ###

# Count data for n labels
count_data <- df_long |>
  na.omit() |>
  group_by(region, condition) |>
  summarise(n = n(), y_pos = max(value, na.rm = TRUE), .groups = "drop")

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/")
pdf("Rg_box_swarm_Sirius.pdf", height = 4.8, width = 6.5)
ggplot(na.omit(df_long), aes(x = region, y = value, color = condition, 
                             group = interaction(condition, region))) +
  geom_boxplot(aes(alpha = NULL, color = NULL), outlier.colour = NA, width = 0.5,
               position = position_dodge(width = 0.75), fill = "white") +
  geom_beeswarm(aes(group = interaction(condition, region)), 
                size = 0.5, dodge.width = 0.75) +
  geom_text(data = count_data, aes(x = region, y = y_pos, label = n,
                                   group = interaction(condition, region), color = NULL),
            position = position_dodge(width = 0.75), size =10, vjust = -0.5,
            show.legend = FALSE, color = "black", alpha = 1) +
  scale_color_manual(values = c("DMSO" = "red", "TX" = "blue")) +
  ylim(0, 0.35) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 2.3, lineend = "round"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line = element_line(colour = "black", linewidth = 2.3, lineend = "round"),
    axis.text = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18)
  ) +
  labs(x = NULL, y = "Value", color = "Condition") + 
  theme(
    axis.text.y = element_text(size = 35, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    
  )

dev.off()




#### Extended Data Fig. 3h - Deff boxplots ###

# Count data for n labels
count_data <- Deff_table_long |>
  na.omit() |>
  group_by(region, condition) |>
  summarise(n = n(), y_pos = max(value, na.rm = TRUE), .groups = "drop")

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_Sirius_FISH/")
pdf("Deff_box_swarm_Sirius.pdf", height = 4.8, width = 5.6)
ggplot(na.omit(Deff_table_long), aes(x = region, y = value, color = condition, 
                                     group = interaction(condition, region))) +
  geom_boxplot(aes(alpha = NULL, color = NULL), outlier.colour = NA, width = 0.5,
               position = position_dodge(width = 0.75), fill = "white") +
  geom_beeswarm(aes(group = interaction(condition, region)), 
                size = 0.5, dodge.width = 0.75) +
  geom_text(data = count_data, aes(x = region, y = y_pos, label = n,
                                   group = interaction(condition, region), color = NULL),
            position = position_dodge(width = 0.75), size = 5, vjust = -0.5,
            show.legend = FALSE, color = "black", alpha = 1) +
  scale_color_manual(values = c("DMSO" = "red", "TX" = "blue")) +
  scale_y_continuous(
    limits = c(0, 0.015),
    breaks = c(0, 0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.014),
    labels = c("0", "", "", "", "", "0.01", "", "")
  ) +
  theme_minimal(base_size = 20) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 1.2, lineend = "round"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line = element_line(colour = "black", linewidth = 1.2, lineend = "round"),
    axis.text = element_text(size = 20, colour = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 18)
  ) +
  labs(x = NULL, y = "Value", color = "Condition") + 
  theme(
    axis.text.y = element_text(size = 23, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    
  )

dev.off()
