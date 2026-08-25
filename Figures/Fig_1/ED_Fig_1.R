library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(txdbmaker)   
library(GenomicRanges)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationHub)


########### Load df_residuals table ########

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


# CPM <- import.bw("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/SLAMseq/coverage/HS37_20250528_50000.bw")
# #bw <- import.bw("file.bw")
CPM <- import.bw("/MyPath/DMSO_20241103_50000.bw")
#

# Convert to data frame (bed-like)
CPM <- as.data.frame(CPM)[, c("seqnames","start","end","score")]
colnames(CPM)[1] <- "chrom"
colnames(CPM)[4] <- "CPM"
CPM$chrom <- paste0("chr", CPM$chrom)
CPM$start <- CPM$start -1



df_residuals <- merge(df_residuals,CPM, by= c("chrom","start","end"), all.x = TRUE)

###### Count genes per bin

genes_gr <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
bins <- fread("/MyPath/hg38_50kb.bed", col.names = c("chrom","start","end"))
bins_gr <- makeGRangesFromDataFrame(bins, starts.in.df.are.0based = TRUE)
bins$gene_count <- countOverlaps(bins_gr, genes_gr)



df_residuals <- merge(df_residuals,bins, by= c("chrom","start","end"), all.x = TRUE)

df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]


######### ED Fig.1a - Violin plot of LOS for 2h_20250328 subset by PC1 sign 

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_range2Mb_K562_2hPD_DpnII_R2_20250326')]
colnames(NT)[6] <- "LOS_residuals"
NT$PC_sign <- ifelse(NT$PC1 > 0, "A","B")

setwd("/MyPath/")

pdf("LOS_2h_20250326_by_PC1_violin.pdf", width=3, height=3)
ggplot(na.omit(NT), aes(x=PC_sign, y=LOS_residuals, fill=PC_sign)) + 
  geom_violin(trim=FALSE, show.legend = FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values = c("red","blue")) +
  ylim(0.4, 0.95) +
  labs(x="Compartment", y = "")+
  geom_hline(yintercept = 0,color="black",linetype=2) +
  guides(fill=guide_legend(title="Compartment"))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=22)) +
  theme(axis.line = element_line(linewidth = 1.4, lineend = "round"),    
        axis.ticks = element_line(linewidth = 1.4, lineend = "round"),  
        axis.ticks.length = unit(0.29, "cm"))
dev.off()


######### ED Fig.1b - Feature heatmap was made using "/sharehome/dlafonta/bin_plot_featureData/bin_normalize_featureData.ipynb"


######### ED Fig.1b - CPM by SPIN state violin



NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','CPM')]

NTTD <- NT[,c(4, 6)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "RIRA"


pdf("CPSp_NT_allSPIN_global_violin.pdf", width=10, height=3)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  scale_y_log10() +
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

######### ED Fig.1d

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_range2Mb_K562_2hPD_DpnII_R2_20250326')]
colnames(NT)[6] <- "LOS_residuals"
NT$PC_sign <- ifelse(NT$PC1 > 0, "A","B")

pdf("LOS_2h_20250326_by_SPIN_alt.pdf", width=7, height=3)
ggplot(na.omit(NT), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(0.4, 0.95) +
  labs(x="SPIN state", y = "Loss Of Structure (LOS)")+
  geom_hline(yintercept = 0,color="black",linetype=2) +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=22)) +
  theme(axis.line = element_line(linewidth = 1.4, lineend = "round"), 
        axis.ticks = element_line(linewidth = 1.4, lineend = "round"),   
        axis.ticks.length = unit(0.29, "cm"))
dev.off()

######### ED Fig. 1e - DpnII-seq violin by SPIN states

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','signal_K562_2hPD_DpnII_R2_20250326')]
colnames(NT)[6] <- "signal"
NT$PC_sign <- ifelse(NT$PC1 > 0, "A","B")
NT <- NT[NT$signal < 5100,]


pdf("DpnIIseq_2h_20250326_by_SPIN_alt.pdf", width=7.45, height=3)
ggplot(na.omit(NT), aes(x=SPIN, y=signal, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(0, 5500) +
  labs(x="SPIN state", y = "")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=22)) +
  theme(axis.line = element_line(linewidth = 1.4, lineend = "round"),  
        axis.ticks = element_line(linewidth = 1.4, lineend = "round"), 
        axis.ticks.length = unit(0.29, "cm"))
dev.off()




######### ED Fig. 1f - Telomere tracks
i <- "chr16"
seg_size <- 20

t <- 1.5
u <- 4.0

x <-df_residuals[df_residuals$chrom == i, ]
pdf("LOS_correction_range2Mb_tracks_20250326_2h_chr16_A_alt.pdf", height = 3.05, width = 9)
par(mfrow=c(2,1), mar=c(0.2, 4, 0, 2) + 0.1,oma=c(5, 0, 0, 0), mgp=c(3, 0.5, 0))
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
axis(1, lwd=u, cex.axis=t, labels=FALSE, tck=-0.075)
axis(2, lwd=u, cex.axis=t, tck=-0.075)
box(bty="l", lwd=u)
dev.off()



t <- 2.3
u <- 4.2

pdf("LOS_correction_range2Mb_tracks_20250326_2h_chr16_B.pdf", height = 7.2, width = 9)
par(mfrow=c(4,1), mar=c(2.7, 4, 1.8, 2) + 0.1,oma=c(5, 0, 0, 0), mgp=c(3, 0.9, 0))
plot(x$start/1000000, x$signal_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(0, 2200), xlab = "", ylab= "", axes=FALSE, lwd=1)
axis(1, lwd=u, cex.axis=t, labels=FALSE) 
axis(2, lwd=u, cex.axis=t)
box(bty="l", lwd=u)


plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326, type="l", col="black", 
     ylim=c(-0.3, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=1)
axis(1, lwd=u, cex.axis=t, labels=FALSE) 
axis(2, lwd=u, cex.axis=t)
abline(h = 0,col = "black", lty=2)
box(bty="l", lwd=u)

plot(x$start/1000000, x$CPM, type="h",
     xlab = "", ylab= "", axes=FALSE,cex.lab=1.2, log = 'y')# ylim=c(0, 100000),)

axis(1, lwd=u, cex.axis=t, labels=FALSE) 
axis(2, lwd=u, cex.axis=t)
     labels = expression(10^-1, 10^3)
box(bty="l", lwd=u)

plot(x$start/1000000, x$gene_count, type="h", col="black", 
     xlab = "", ylab= "", axes=FALSE, lwd=1)#,ylim=c(-0.2, 0.1))
axis(1, lwd=u, cex.axis=t, labels=TRUE,mgp=c(3, 1.4, 0)) 
axis(2, lwd=u, cex.axis=t)
box(bty="l", lwd=u)

dev.off()




##############################################################################################################
######################   Telomere analysis  ##################################################################


###### Subset telomeres 

chrom_sizes <- df_residuals %>%
  group_by(chrom) %>%
  summarise(size = max(end)) %>%
  ungroup()

size <- chrom_sizes$size


telomere_bins_p <- df_residuals %>%
  left_join(chrom_sizes, by = "chrom") %>%
  filter(
    start <= 10000000) %>%
  dplyr::select(-size)

telomere_bins_q <- df_residuals %>%
  left_join(chrom_sizes, by = "chrom") %>%
  filter(end >= (size - 10000000)  # Within 10 Mb of q-arm telomere
  ) %>%
  dplyr::select(-size)


# Count occurrences of "Speckle" in SPIN column grouped by chrom
result_p <- telomere_bins_p %>%
  group_by(chrom) %>%
  summarise(
    speckle_bin_p = sum(str_count(SPIN, "Speckle"), na.rm = TRUE)
  )

result_q <- telomere_bins_q %>%
  group_by(chrom) %>%
  summarise(
    speckle_bin_q = sum(str_count(SPIN, "Speckle"), na.rm = TRUE)
  )

result <- merge(result_p,result_q, by="chrom")




telomere_bins <- df_residuals %>%
  left_join(chrom_sizes, by = "chrom") %>%
  filter(
    start <= 10000000 |  # Within 10 Mb of p-arm telomere
      end >= (size - 10000000)  # Within 10 Mb of q-arm telomere
  ) %>%
  dplyr::select(-size)

# Calculate counts in telomere bins and genome-wide
# Remove NAs from continuous column
telomere_bins <- telomere_bins %>%
  filter(!is.na(LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326))
telomere_bins <- telomere_bins %>%
  filter(!is.na(SPIN))

bg <- df_residuals %>%
  filter(!is.na(LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326))
bg <- df_residuals %>%
  filter(!is.na(SPIN))

telomere_counts <- telomere_bins %>%
  count(SPIN, name = "telomere_n")

genome_counts <- bg %>%
  count(SPIN, name = "genome_n")

SPINbin_counts_telomere <- telomere_counts %>%
  left_join(genome_counts, by = "SPIN") %>%
  mutate(
    n = telomere_n,
    percent = round(telomere_n / genome_n * 100, 1)
  )

########


######### ED Fig. 1g - Telomere violins

pdf("LOSres_2h_20250326_by_SPIN_atTelomeres_alt.pdf", width=7, height=3)
ggplot(telomere_bins, aes(x = SPIN, y = LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326, 
                          fill = SPIN)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white")+
  geom_text(data = SPINbin_counts_telomere, 
            aes(x = SPIN, y = Inf, 
                label = paste0("n=", n, "\n", percent, "%")),
            vjust = 1.2, inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.36, 0.26) +
  labs(x="SPIN state", y = "Loss Of Structure (LOS)")+
  geom_hline(yintercept = 0,color="black",linetype=2) +
  theme_classic()+
  theme(axis.text=element_text(size=25))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.line = element_line(linewidth = 1.4, lineend = "round"),   
        axis.ticks = element_line(linewidth = 1.4, lineend = "round"),
        axis.ticks.length = unit(0.15, "cm"))


dev.off()



############################# Centromere analysis  ###############################



ah <- AnnotationHub()

# Search for centromere data
query(ah, c("centromere", "Homo sapiens"))

# Get centromere positions for hg38
centromeres <- ah[["AH107354"]]

# Convert to data frame
centromeres_single <- range(split(centromeres, seqnames(centromeres)))
centromeres_single <- unlist(centromeres_single)


flank_size <- 10000000  

peri_centromeric <- as.data.frame(centromeres_single)
peri_centromeric$start  <- peri_centromeric$start - flank_size
peri_centromeric$end  <- peri_centromeric$end + flank_size

rownames(peri_centromeric) <- NULL
colnames(peri_centromeric)[1] <- "chrom"

peri_centromeric_gr <- makeGRangesFromDataFrame(peri_centromeric, keep.extra.columns = TRUE)
df_residuals_gr <- makeGRangesFromDataFrame(df_residuals, keep.extra.columns = TRUE)
peri_centromeric_subset <- subsetByOverlaps(df_residuals_gr, peri_centromeric_gr)
peri_centromeric_subset_df <- as.data.frame(peri_centromeric_subset)



# Calculate counts in telomere bins and genome-wide
# Remove NAs from continuous column
peri_centromeric_subset_df <- peri_centromeric_subset_df %>%
  filter(!is.na(LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326))
peri_centromeric_subset_df <- peri_centromeric_subset_df %>%
  filter(!is.na(SPIN))

bg <- df_residuals %>%
  filter(!is.na(LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326))
bg <- df_residuals %>%
  filter(!is.na(SPIN))

centromere_counts <- peri_centromeric_subset_df %>%
  count(SPIN, name = "centromere_n")

genome_counts <- bg %>%
  count(SPIN, name = "genome_n")

SPINbin_counts_centromere <- centromere_counts %>%
  left_join(genome_counts, by = "SPIN") %>%
  mutate(
    n = centromere_n,
    percent = round(centromere_n / genome_n * 100, 1)
  )



######### ED Fig. 1h - Centromere Violin

#Violin plotting the LOS residuals at centromeres

pdf("LOSres_2h_20250326_by_SPIN_atCentromeres_alt.pdf", width=7, height=3)
ggplot(peri_centromeric_subset_df, aes(x = SPIN, y = LOS_residuals_range2Mb_K562_2hPD_DpnII_R2_20250326, 
                          fill = SPIN)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white")+
  geom_text(data = SPINbin_counts_centromere, 
            aes(x = SPIN, y = Inf, 
                label = paste0("n=", n, "\n", percent, "%")),
            vjust = 1.2, inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.36, 0.26) +
  labs(x="SPIN state", y = "Loss Of Structure (LOS)")+
  geom_hline(yintercept = 0,color="black",linetype=2) +
  theme_classic()+
  theme(axis.text=element_text(size=25))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.line = element_line(linewidth = 1.4, lineend = "round"), 
        axis.ticks = element_line(linewidth = 1.4, lineend = "round"), 
        axis.ticks.length = unit(0.15, "cm"))

dev.off()



######## ED Fig. 2h  -   Growth curves for long culture without media addition  ###########
growth_data <- read.table("/MyPath/20240321_K562_long culture_R1.tsv",
                          sep="\t", header=TRUE)

###### plot growth curves  ####

growth_data$M_cellsPerMl <- growth_data$AVG.cells.ml/1000000


ylim.prim <- c(0, 1.3)  
ylim.sec <- c(0, 100)  

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]



pdf("Growth_curves_longCulture.pdf", height = 3, width = 7.1)
ggplot(growth_data, aes(Days, M_cellsPerMl)) +
  geom_point() +
  geom_smooth(
    method = "lm",
    formula = y ~ splines::ns(x, df = 9), 
    se = FALSE 
  ) +
  geom_smooth(aes(y = a + AVG..live*b),method = "lm",
              formula = y ~ splines::ns(x, df = 9),
              se = FALSE 
              ,color = "red") +
  geom_point(aes(y = a + AVG..live*b) ,color = "red") +
  theme_classic()+
  scale_y_continuous("", sec.axis = sec_axis(~ (. - a)/b, name = "")) +
  scale_x_continuous("", breaks = seq(0,21, by=5)) +
  theme(axis.title.y.right = element_text(color = "red")) +
  theme(axis.title.y.left = element_text(color = "blue")) +
  theme(axis.text=element_text(size=23))+
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.line = element_line(linewidth = 1.3, lineend = "round"),
        axis.ticks = element_line(linewidth = 1.3, lineend = "round"),
        axis.ticks.length = unit(0.15, "cm")) +
  theme(axis.ticks.length = unit(0.15, "cm"),
        axis.text.x = element_text(margin = margin(t = 5)))
dev.off()



## ED Fig. 1k - Plot track of fresh, old and HB R1

seg_size <- 25

begin <- 55000000
term <- 85000000

t <- 1.3
u <- 2.75

x <- df_residuals[df_residuals$chrom == "chr11",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_Growth_track_ch11_wideA_alt2.pdf", height = 2.7, width = 4.7)
par(mfrow=c(2,1), mar=c(0.2, 4, 0, 2) + 0.1,oma=c(5, 0, 0, 0), mgp=c(3, 0.5, 0))

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
     xlab = "", ylab= "", ylim=c(-1, 1), axes=FALSE)
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
axis(1, lwd=u, cex.axis=t, labels=FALSE, tck=-0.075)
axis(2, lwd=u, cex.axis=t, tck=-0.075)
box(bty="l", lwd=u)

dev.off()


t <- 2.3
u <- 3.2

x <- df_residuals[df_residuals$chrom == "chr11",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_Growth_track_ch11_wideB_alt2.pdf", height = 5.1, width = 4.7)
par(mfrow=c(3,1), mar=c(2.7, 4, 1.8, 2) + 0.1,oma=c(5, 0, 0, 0), mgp=c(3, 1.8, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_contCD_120mPD_DpnII_20240327, type="l", col="black",
     ylim=c(-0.425, 0.225), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=u, cex.axis=t, labels=FALSE)
axis(2, lwd=u, cex.axis=t, tck=-0.08, mgp=c(3, 1, 0))
box(bty="l", lwd=u)

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_old120mPD_DpnII_R2_20240327, type="l", col="black",
     ylim=c(-0.425, 0.225), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=u, cex.axis=t, labels=FALSE)
axis(2, lwd=u, cex.axis=t, tck=-0.08, mgp=c(3, 1, 0))
box(bty="l", lwd=u)

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_120mPD_DpnII_R1_20170116, type="l", col="black",
     ylim=c(-0.425, 0.225), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=u, cex.axis=2, labels=TRUE, mgp=c(3, 1.3, 0))
axis(2, lwd=u, cex.axis=t, tck=-0.08, mgp=c(3, 1, 0))
box(bty="l", lwd=u)

dev.off()



## ED Fig. 1k - Plot violins for fresh, old and HB R1 LOS residuals

################  Global violins    ########################

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562_contCD_120mPD_DpnII_20240327','LOS_residuals_range2Mb_K562_old120mPD_DpnII_R2_20240327','LOS_residuals_range2Mb_K562_120mPD_DpnII_R1_20170116')]


NTTD <- NT[,c(4, 6)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "fresh"
NTTD <- NTTD[NTTD$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN_global_fresh_violin_alt2.pdf", width=4.5, height=2.05)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61","#ABDDA4",  "#5E4FA2")) +
  ylim(-0.42, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2.2, lineend = "round"),    
        axis.ticks = element_line(linewidth = 2.2, lineend = "round"),  
        axis.ticks.length = unit(0.29, "cm"))

dev.off()


U14h <- NT[,c(4, 7)]
colnames(U14h)[2] <- "LOS_residuals"
U14h$treat <- "old"
U14h <- U14h[U14h$SPIN == c("Speckle", "Interior_Act1", "Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN_global_old_violin_alt2.pdf", width=4.5, height=2.05)
ggplot(na.omit(U14h), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43","#FDAE61", "#ABDDA4",  "#5E4FA2")) +
  ylim(-0.42, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2.2, lineend = "round"),  
        axis.ticks = element_line(linewidth = 2.2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()



RIRA <- NT[,c(4, 8)]
colnames(RIRA)[2] <- "LOS_residuals"
RIRA$treat <- "RNase A"
RIRA <- RIRA[RIRA$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN_global_HBR1_violin_alt2.pdf", width=4.5, height=2.05)
ggplot(na.omit(RIRA), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61","#ABDDA4",  "#5E4FA2")) +
  ylim(-0.42, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2.2, lineend = "round"), 
        axis.ticks = element_line(linewidth = 2.2, lineend = "round"),    
        axis.ticks.length = unit(0.29, "cm"))

dev.off()



