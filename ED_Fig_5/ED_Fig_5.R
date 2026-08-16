library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)
library(clusterProfiler)


df_residuals <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/dataframes/df_residuals_20260812_50kb.bed",
                           sep="\t", header=TRUE)
df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]

df_residuals <- df_residuals %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina", "Lamina_Like"),
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
df_residuals$Lamina_Like_seg <- ifelse(df_residuals$SPIN == "Lamina_Like",df_residuals$Lamina_Like_seg, 10)

df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]


df_residuals_250kb <-  read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/250kb/dataframes/df_residuals_20260225_250kb.bed",
                                  sep="\t", header=TRUE)
df_residuals_250kb <- df_residuals_250kb[order(df_residuals_250kb$chrom, df_residuals_250kb$start,df_residuals_250kb$end),]
df_residuals_250kb <- df_residuals_250kb %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina", "Lamina_Like"),
                  ordered = TRUE))

###############################  GO analysis for exonic reads  #######################################################################################

###################### GO analysis setting p-adjusted to 0.05  #########################


DGE_HS42_V1_allTCReads <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/SLAMseq/3_5_26/NextSeq_6_12_25_DGE_all_counts_allgenes.tsv",
                                     sep="\t", header=TRUE)
DGE_HS42_V1_allTCReads$treatment <- replace(DGE_HS42_V1_allTCReads$treatment, DGE_HS42_V1_allTCReads$treatment == "HS42vscontrol", "Heat_shock_80min_allReads")
DGE_HS42_V1_allTCReads$treatment <- replace(DGE_HS42_V1_allTCReads$treatment, DGE_HS42_V1_allTCReads$treatment == "V1vsDMSO", "Speckle-dTAGvsDMSO_6h_allReads")
DGE_HS42_V1_allTCReads$treatment <- replace(DGE_HS42_V1_allTCReads$treatment, DGE_HS42_V1_allTCReads$treatment == "V1vsNeg", "Speckle-dTAGvsNEG_6h_allReads")

DGE_HS42_V1_exonTCReads <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/SLAMseq/3_5_26/NextSeq_6_12_25_DGE_exonic_counts_allgenes.tsv",
                                      sep="\t", header=TRUE)
DGE_HS42_V1_exonTCReads$treatment <- replace(DGE_HS42_V1_exonTCReads$treatment, DGE_HS42_V1_exonTCReads$treatment == "HS42vscontrol", "Heat_shock_80min_exonsOnly")
DGE_HS42_V1_exonTCReads$treatment <- replace(DGE_HS42_V1_exonTCReads$treatment, DGE_HS42_V1_exonTCReads$treatment == "V1vsDMSO", "Speckle-dTAGvsDMSO_6h_exonsOnly")
DGE_HS42_V1_exonTCReads$treatment <- replace(DGE_HS42_V1_exonTCReads$treatment, DGE_HS42_V1_exonTCReads$treatment == "V1vsNeg", "Speckle-dTAGvsNEG_6h_exonsOnly")

DGE_HS42_V1_intronTCReads <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/SLAMseq/3_5_26/NextSeq_6_12_25_DGE_intronic_counts_allgenes.tsv",
                                        sep="\t", header=TRUE)
DGE_HS42_V1_intronTCReads$treatment <- replace(DGE_HS42_V1_intronTCReads$treatment, DGE_HS42_V1_intronTCReads$treatment == "HS42vscontrol", "Heat_shock_80min_intronssOnly")
DGE_HS42_V1_intronTCReads$treatment <- replace(DGE_HS42_V1_intronTCReads$treatment, DGE_HS42_V1_intronTCReads$treatment == "V1vsDMSO", "Speckle-dTAGvsDMSO_6h_intronssOnly")
DGE_HS42_V1_intronTCReads$treatment <- replace(DGE_HS42_V1_intronTCReads$treatment, DGE_HS42_V1_intronTCReads$treatment == "V1vsNeg", "Speckle-dTAGvsNEG_6h_intronssOnly")

DGE_allTCReads <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/SLAMseq/3_5_26/NextSeq_12_30_24_DGE_all_counts_allgenes.tsv",
                             sep="\t", header=TRUE)
DGE_allTCReads$treatment <- replace(DGE_allTCReads$treatment, DGE_allTCReads$treatment == "treatment_triptolide_DRB_vs_control_triptolide_DRB", "Transcription_block_allReads")
DGE_allTCReads$treatment <- replace(DGE_allTCReads$treatment, DGE_allTCReads$treatment == "treatment_AMO_4h_30uM_vs_control_AMO_4h_30uM", "U1AMO_4h_30uM_allReads")
DGE_allTCReads$treatment <- replace(DGE_allTCReads$treatment, DGE_allTCReads$treatment == "treatment_AMO_4h_60uM_vs_control_AMO_4h_60uM", "U1AMO_4h_60uM_allReads")
DGE_allTCReads$treatment <- replace(DGE_allTCReads$treatment, DGE_allTCReads$treatment == "treatment_AMO_8h_30uM_vs_control_AMO_8h_30uM", "U1AMO_8h_30uM_allReads")
DGE_allTCReads$treatment <- replace(DGE_allTCReads$treatment, DGE_allTCReads$treatment == "treatment_AMO_8h_60uM_vs_control_AMO_8h_60uM", "U1AMO_8h_60uM_allReads")

DGE_exonTCReads <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/SLAMseq/3_5_26/NextSeq_12_30_24_DGE_exonic_counts_allgenes.tsv",
                              sep="\t", header=TRUE)
DGE_exonTCReads$treatment <- replace(DGE_exonTCReads$treatment, DGE_exonTCReads$treatment == "treatment_triptolide_DRB_vs_control_triptolide_DRB", "Transcription_block_exonsOnly")
DGE_exonTCReads$treatment <- replace(DGE_exonTCReads$treatment, DGE_exonTCReads$treatment == "treatment_AMO_4h_30uM_vs_control_AMO_4h_30uM", "U1AMO_4h_30uM_exonsOnly")
DGE_exonTCReads$treatment <- replace(DGE_exonTCReads$treatment, DGE_exonTCReads$treatment == "treatment_AMO_4h_60uM_vs_control_AMO_4h_60uM", "U1AMO_4h_60uM_exonsOnly")
DGE_exonTCReads$treatment <- replace(DGE_exonTCReads$treatment, DGE_exonTCReads$treatment == "treatment_AMO_8h_30uM_vs_control_AMO_8h_30uM", "U1AMO_8h_30uM_exonsOnly")
DGE_exonTCReads$treatment <- replace(DGE_exonTCReads$treatment, DGE_exonTCReads$treatment == "treatment_AMO_8h_60uM_vs_control_AMO_8h_60uM", "U1AMO_8h_60uM_exonsOnly")

DGE_intronTCReads <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/SLAMseq/3_5_26/NextSeq_12_30_24_DGE_intronic_counts_allgenes.tsv",
                                sep="\t", header=TRUE)
DGE_intronTCReads$treatment <- replace(DGE_intronTCReads$treatment, DGE_intronTCReads$treatment == "treatment_triptolide_DRB_vs_control_triptolide_DRB", "Transcription_block_intronsOnly")
DGE_intronTCReads$treatment <- replace(DGE_intronTCReads$treatment, DGE_intronTCReads$treatment == "treatment_AMO_4h_30uM_vs_control_AMO_4h_30uM", "U1AMO_4h_30uM_intronsOnly")
DGE_intronTCReads$treatment <- replace(DGE_intronTCReads$treatment, DGE_intronTCReads$treatment == "treatment_AMO_4h_60uM_vs_control_AMO_4h_60uM", "U1AMO_4h_60uM_intronsOnly")
DGE_intronTCReads$treatment <- replace(DGE_intronTCReads$treatment, DGE_intronTCReads$treatment == "treatment_AMO_8h_30uM_vs_control_AMO_8h_30uM", "U1AMO_8h_30uM_intronsOnly")
DGE_intronTCReads$treatment <- replace(DGE_intronTCReads$treatment, DGE_intronTCReads$treatment == "treatment_AMO_8h_60uM_vs_control_AMO_8h_60uM", "U1AMO_8h_60uM_intronsOnly")


##### Consolidate tables ###
### All
cols_keep <- c("Gene", "FoldChange", "log2FoldChange", "padj", "Color_factor", "treatment")

DGE_HS42_V1_allTCReads <- DGE_HS42_V1_allTCReads[, cols_keep]
DGE_allTCReads         <- DGE_allTCReads[, cols_keep]

DGE_allTCReads_combined <- rbind(DGE_HS42_V1_allTCReads, DGE_allTCReads)

### Exons
DGE_HS42_V1_exonTCReads <- DGE_HS42_V1_exonTCReads[, cols_keep]
DGE_exonTCReads         <- DGE_exonTCReads[, cols_keep]

DGE_exonTCReads_combined <- rbind(DGE_HS42_V1_exonTCReads, DGE_exonTCReads)

### Introns
DGE_HS42_V1_intronTCReads <- DGE_HS42_V1_intronTCReads[, cols_keep]
DGE_intronTCReads         <- DGE_intronTCReads[, cols_keep]

DGE_intronTCReads_combined <- rbind(DGE_HS42_V1_intronTCReads, DGE_intronTCReads)

##### Make Log2FC column
DGE_allTCReads_combined$Log2FC <- log2(DGE_allTCReads_combined$FoldChange)
DGE_exonTCReads_combined$Log2FC <- log2(DGE_exonTCReads_combined$FoldChange)
DGE_intronTCReads_combined$Log2FC <- log2(DGE_intronTCReads_combined$FoldChange)


#### Make significance column

DGE_allTCReads_combined$Significance <- ifelse(DGE_allTCReads_combined$Color_factor %in% c("Significantly_Up", "Significantly_Down"), 
                                               "Significant", 
                                               "Not_Significant")

DGE_exonTCReads_combined$Significance <- ifelse(DGE_exonTCReads_combined$Color_factor %in% c("Significantly_Up", "Significantly_Down"), 
                                                "Significant", 
                                                "Not_Significant")

DGE_intronTCReads_combined$Significance <- ifelse(DGE_intronTCReads_combined$Color_factor %in% c("Significantly_Up", "Significantly_Down"), 
                                                  "Significant", 
                                                  "Not_Significant")


######## Extended Data Fig 4c



###################   Plotting Split Violins  ###############

group_specs <- list(
  list(data = DGE_allTCReads_combined,    treatment_src = "Heat_shock_80min_allReads",    label = "All Reads"),
  list(data = DGE_exonTCReads_combined, treatment_src = "Heat_shock_80min_exonsOnly", label = "Exons"),
  list(data = DGE_intronTCReads_combined, treatment_src = "Heat_shock_80min_intronssOnly", label = "Introns")
)

# Build a unified long table with a new "group" column
combined_all <- bind_rows(lapply(group_specs, function(gs) {
  gs$data %>%
    filter(treatment == gs$treatment_src) %>%
    mutate(group = gs$label)
}))

group_levels <- c("All Reads", "Exons", "Introns") 
combined_all$group <- factor(combined_all$group, levels = group_levels)

scale_factor <- 0.2
box_width    <- 0.04

compute_half_violin <- function(df, side, group_col = "group") {
  groups <- levels(df[[group_col]])
  bind_rows(lapply(groups, function(grp) {
    vals <- df$Log2FC[df[[group_col]] == grp]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 2) return(NULL)
    d         <- density(vals)
    dens_norm <- d$y / max(d$y)
    dens      <- if (side == "left") -dens_norm else dens_norm
    data.frame(y = d$x, density = dens, group = grp, side = side)
  }))
}

compute_boxplot_stats <- function(df, side, group_col = "group") {
  groups <- levels(df[[group_col]])
  bind_rows(lapply(groups, function(grp) {
    vals     <- df$Log2FC[df[[group_col]] == grp]
    vals     <- vals[!is.na(vals)]
    x_center <- as.numeric(factor(grp, levels = group_levels))
    x_offset <- if (side == "left") -scale_factor * 0.25 else scale_factor * 0.25
    data.frame(
      group  = grp,
      x      = x_center + x_offset,
      lower  = quantile(vals, 0.25),
      middle = median(vals),
      upper  = quantile(vals, 0.75),
      ymin   = min(vals[vals >= quantile(vals, 0.25) - 1.5 * IQR(vals)]),
      ymax   = max(vals[vals <= quantile(vals, 0.75) + 1.5 * IQR(vals)])
    )
  }))
}

compute_annotations <- function(df, group_col = "group") {
  groups <- levels(df[[group_col]])
  bind_rows(lapply(groups, function(grp) {
    sub      <- df[df[[group_col]] == grp & !is.na(df$Log2FC), ]
    n_total  <- nrow(sub)
    n_up     <- sum(sub$Color_factor == "Significantly_Up",   na.rm = TRUE)
    n_down   <- sum(sub$Color_factor == "Significantly_Down", na.rm = TRUE)
    data.frame(
      group = grp,
      x     = as.numeric(factor(grp, levels = group_levels)),
      y     = 3,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  }))
}


sig_all <- combined_all %>% filter(Significance == "Significant")

left_dens  <- compute_half_violin(combined_all, "left")
right_dens <- compute_half_violin(sig_all,      "right")

left_dens$x_pos  <- as.numeric(factor(left_dens$group,  levels = group_levels)) + left_dens$density  * scale_factor
right_dens$x_pos <- as.numeric(factor(right_dens$group, levels = group_levels)) + right_dens$density * scale_factor

box_left  <- compute_boxplot_stats(combined_all, "left")
box_right <- compute_boxplot_stats(sig_all,      "right")

annot_df  <- compute_annotations(combined_all)
annot_df  <- compute_annotations(combined_all)


setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")

pdf("HS_SplitViolin_SLAMseq_combined.pdf", height = 5, width = 9)

ggplot() +
  geom_polygon(data = left_dens,
               aes(x = x_pos, y = y, group = group, fill = "All"),
               alpha = 0.7, color = "black", linewidth = 0.8) +
  geom_polygon(data = right_dens,
               aes(x = x_pos, y = y, group = group, fill = "Significant"),
               alpha = 0.7, color = "black", linewidth = 0.8) +
  geom_boxplot(data = box_left,
               aes(x = x, ymin = ymin, lower = lower, middle = middle,
                   upper = upper, ymax = ymax, group = group),
               stat = "identity", width = box_width, fill = "white", alpha = 0.8) +
  geom_boxplot(data = box_right,
               aes(x = x, ymin = ymin, lower = lower, middle = middle,
                   upper = upper, ymax = ymax, group = group),
               stat = "identity", width = box_width, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = annot_df, aes(x = x, y = 4.9, label = label),
            vjust = -0.3, size = 8, lineheight = 0.9) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    limits = c(0.7, length(group_levels) + 0.3)
  ) +
  scale_fill_manual(values = c("All" = "#7fbfff", "Significant" = "#ff7f7f")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 24, margin = margin(t = 5)),
    axis.title.x     = element_blank(),
    axis.ticks.x     = element_line(linewidth = 2.7, lineend = "round"),
    axis.ticks.y     = element_line(linewidth = 2.7, lineend = "round"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line        = element_line(linewidth = 2.7, lineend = "round"),
    axis.text.y      = element_text(size = 37, margin = margin(r = 5)),
    axis.title.y     = element_blank(),
    legend.position  = "none"
  ) +
  labs(y = "Log2 Fold Change")

dev.off()



######## Extended Data Figure 5d

###############################  GO analysis for exonic reads  #######################################################################################

###################### GO analysis setting p-adjusted to 0.05  #########################
m <- useEnsembl(biomart = "genes",
                dataset = "hsapiens_gene_ensembl",
                mirror  = "useast")

symb <- keys(org.Hs.eg.db, "SYMBOL")
attrs  <- listAttributes(mart = m)
df <- getBM(c("ensembl_gene_id","hgnc_symbol"), "hgnc_symbol", symb, m)
colnames(df)[2] <- 'geneName' 

colnames(DGE_allTCReads_combined)[colnames(DGE_allTCReads_combined) == "Gene"] <- "ensembl_gene_id"
DGE_allTCReads_combined <- merge(DGE_allTCReads_combined,df, by="ensembl_gene_id", all.x= TRUE) 
up_genes_SLAM <- DGE_allTCReads_combined[(DGE_allTCReads_combined$Color_factor == "Significantly_Up") & DGE_allTCReads_combined$treatment == "Heat_shock_80min_allReads",]

up_genes <- up_genes_SLAM[['ensembl_gene_id']]



ego_SLAM <- enrichGO(gene= up_genes,
                     OrgDb= org.Hs.eg.db,
                     keyType = 'ENSEMBL',
                     ont= "all")

s_ego_SLAM <- clusterProfiler::simplify(ego_SLAM)

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/GO_DotPlot")
pdf("GOenrichUP_SLAMseq_DotPlot_exon_simple.pdf", height = 6, width = 7)
dotplot(s_ego_SLAM, showCategory = 7, font.size = 12, label_format = 25) +
  theme_minimal() +
  theme(
    axis.text.y  = element_text(size = 16), 
    axis.text.x  = element_text(size = 16), 
    axis.title   = element_text(size = 16), 
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 12),
    plot.title   = element_text(size = 15, face = "bold")
  )
dev.off()




######## Ext Data Fig 5f



setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/tracks")
seg_size <- 25

begin <- 33000000
term <- 65000000

t <- 1.3
lw  <- 2.6


x <- df_residuals[df_residuals$chrom == "chr3",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_HS_CdSulf_track_ch3_zoom_A_alt.pdf", height = 2.9, width = 4)
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
axis(1, lwd=lw , cex.axis=t, labels=FALSE, tck=-0.075)
axis(2, lwd=lw , cex.axis=t, tck=-0.075)
box(bty="l", lwd=lw )

dev.off()




t <- 1.25
lw  <- 2.8
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/tracks")

x <- df_residuals[df_residuals$chrom == "chr3",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_HSonly_track_ch3_zoom_B.pdf", height = 3, width = 4)
par(mfrow=c(2,1), mar=c(0.75, 4, 1.8, 2) + 0.1,oma=c(1, 0, 0, 0), mgp=c(3, 0.6, 0))
plot(x$start/1000000, x$LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif, type="l", col="black",
     ylim=c(-0.22, 0.15), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=t, tck=-0.08,labels=TRUE)
axis(2, lwd=lw , cex.axis=t, tck=-0.08, mgp=c(3, 0.6, 0))
box(bty="l", lwd=lw )


dev.off()



## Extended Data Fig.5f - Global violins 

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/violin")

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif','LOS_residuals_range2Mb_K562_CdSulfate_120mPD_DpnII_20240327_dif')]


HS42 <- NT[,c(4, 6)]
colnames(HS42)[2] <- "LOS_residuals"
HS42$treat <- "TPL/DRB"
HS42 <- HS42[HS42$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN_global_HS42_violin.pdf", width=4, height=2.05)
ggplot(na.omit(HS42), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#ABDDA4", "#5E4FA2")) +
  ylim(-0.215, 0.145) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2, lineend = "round"), 
        axis.ticks = element_line(linewidth = 2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

options(scipen = 999)
aggregate(LOS_residuals ~ SPIN, data = na.omit(HS42), FUN = mean)









