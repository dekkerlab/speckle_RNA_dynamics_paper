library(tidyverse)
library(httr)
library(jsonlite)
library(GenomicRanges)  

df_residuals <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/dataframes/df_residuals_20260812_50kb.bed",
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



df_residuals_150 <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/150kb/dataframes/df_residuals_20260604_150kb.bed",
                           sep="\t", header=TRUE)
df_residuals_150$chrom <- factor(df_residuals_150$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals_150 <- df_residuals_150[order(df_residuals_150$chrom, df_residuals_150$start,df_residuals_150$end),]

df_residuals_150 <- df_residuals_150 %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),
                  ordered = TRUE))

df_residuals_150$seg_1 <- -1
df_residuals_150 <- cbind(df_residuals_150, replicate(8, df_residuals_150$seg_1))
n <- ncol(df_residuals_150)
names(df_residuals_150)[(n-7):n] <- c("Speckle_seg", "Interior_Act1_seg", "Interior_Act2_seg",  "Interior_Repr1_seg","Interior_Repr2_seg","Near_Lm1_seg","Near_Lm2_seg","Lamina_seg")


df_residuals_150$Speckle_seg <- ifelse(df_residuals_150$SPIN == "Speckle",df_residuals_150$Speckle_seg, 10)
df_residuals_150$Interior_Act1_seg <- ifelse(df_residuals_150$SPIN == "Interior_Act1",df_residuals_150$Interior_Act1_seg, 10)
df_residuals_150$Interior_Act2_seg <- ifelse(df_residuals_150$SPIN == "Interior_Act2",df_residuals_150$Interior_Act2_seg, 10)
df_residuals_150$Interior_Repr1_seg <- ifelse(df_residuals_150$SPIN == "Interior_Repr1",df_residuals_150$Interior_Repr1_seg, 10)
df_residuals_150$Interior_Repr2_seg <- ifelse(df_residuals_150$SPIN == "Interior_Repr2",df_residuals_150$Interior_Repr2_seg, 10)
df_residuals_150$Near_Lm1_seg <- ifelse(df_residuals_150$SPIN == "Near_Lm1",df_residuals_150$Near_Lm1_seg, 10)
df_residuals_150$Near_Lm2_seg <- ifelse(df_residuals_150$SPIN == "Near_Lm2",df_residuals_150$Near_Lm2_seg, 10)
df_residuals_150$Lamina_seg <- ifelse(df_residuals_150$SPIN == "Lamina",df_residuals_150$Lamina_seg, 10)

df_residuals_150$chrom <- factor(df_residuals_150$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals_150 <- df_residuals_150[order(df_residuals_150$chrom, df_residuals_150$start,df_residuals_150$end),]



############### Upload SLAM-seq data 

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




###################   Plotting Split Violins  ###############
#### Fig. 4c


group_specs <- list(
  list(data = DGE_allTCReads_combined,    treatment_src = "Speckle-dTAGvsNEG_6h_allReads",    label = "All Reads"),
  list(data = DGE_intronTCReads_combined, treatment_src = "Speckle-dTAGvsNEG_6h_intronssOnly", label = "Introns")
)

# Build a unified long table with a new "group" column
combined_all <- bind_rows(lapply(group_specs, function(gs) {
  gs$data %>%
    filter(treatment == gs$treatment_src) %>%
    mutate(group = gs$label)
}))

group_levels <- c("All Reads", "Introns")   # x-axis order
combined_all$group <- factor(combined_all$group, levels = group_levels)

sig_all <- combined_all %>% filter(Significance == "Significant")



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
      y     = 4,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  }))
}



left_dens  <- compute_half_violin(combined_all, "left")
right_dens <- compute_half_violin(sig_all,      "right")

left_dens$x_pos  <- as.numeric(factor(left_dens$group,  levels = group_levels)) + left_dens$density  * scale_factor
right_dens$x_pos <- as.numeric(factor(right_dens$group, levels = group_levels)) + right_dens$density * scale_factor

box_left  <- compute_boxplot_stats(combined_all, "left")
box_right <- compute_boxplot_stats(sig_all,      "right")

annot_df  <- compute_annotations(combined_all)


setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")

pdf("SpdTAGvsNeg_SplitViolin_SLAMseq_combined.pdf", height = 5, width = 7)

lw <- 1.7
ft <- 32

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
  geom_text(data = annot_df, aes(x = x, y = 7, label = label),
            vjust = -0.3, size = 8, lineheight = 0.9) +
  scale_x_continuous(
    breaks = seq_along(group_levels),
    labels = group_levels,
    limits = c(0.7, length(group_levels) + 0.3)
  ) +
  scale_fill_manual(values = c("All" = "#7fbfff", "Significant" = "#ff7f7f")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(ylim = c(-7, 9.5)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x      = element_text(size = 24, margin = margin(t = 5)),
    axis.title.x     = element_blank(),
    axis.ticks.x     = element_line(linewidth = lw, lineend = "round"),
    axis.ticks.y     = element_line(linewidth = lw, lineend = "round"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line        = element_line(linewidth = lw, lineend = "round"),
    axis.text.y      = element_text(size = ft, margin = margin(r = 5)),
    axis.title.y     = element_blank(),
    legend.position  = "none"
  ) +
  labs(y = "Log2 Fold Change")

dev.off()



###########  Extended Data Fig. 4a

# 
# ## Helper functions 
compute_half_violin <- function(df, side, group_col = "treatment") {
  treatments <- unique(df[[group_col]])
  result <- lapply(treatments, function(trt) {
    vals <- df$Log2FC[df[[group_col]] == trt]
    vals <- vals[!is.na(vals)]
    d <- density(vals)
    dens_norm <- d$y / max(d$y)
    dens <- if (side == "left") -dens_norm else dens_norm
    data.frame(y = d$x, density = dens, treatment = trt, side = side)
  })
  bind_rows(result)
}
# 
# # Update x_offset in compute_boxplot_stats to match scale_factor
compute_boxplot_stats <- function(df, side, group_col = "treatment") {
  treatments <- unique(df[[group_col]])
  result <- lapply(treatments, function(trt) {
    vals <- df$Log2FC[df[[group_col]] == trt]
    vals <- vals[!is.na(vals)]
    x_center <- as.numeric(factor(trt, levels = treatment_levels))
    x_offset <- if (side == "left") -scale_factor * 0.25 else scale_factor * 0.25  
    data.frame(
      treatment = trt,
      x      = x_center + x_offset,
      lower  = quantile(vals, 0.25),
      middle = median(vals),
      upper  = quantile(vals, 0.75),
      ymin   = min(vals[vals >= quantile(vals, 0.25) - 1.5 * IQR(vals)]),
      ymax   = max(vals[vals <= quantile(vals, 0.75) + 1.5 * IQR(vals)])
    )
  })
  bind_rows(result)
}
# 
compute_annotations <- function(all_df, group_col = "treatment") {
  treatments <- unique(all_df[[group_col]])
  result <- lapply(treatments, function(trt) {
    all_vals <- all_df$Log2FC[all_df[[group_col]] == trt & !is.na(all_df$Log2FC)]
    n_total  <- length(all_vals)
    n_up     <- sum(all_df[[group_col]] == trt & all_df$Color_factor == "Significantly_Up",   na.rm = TRUE)
    n_down   <- sum(all_df[[group_col]] == trt & all_df$Color_factor == "Significantly_Down", na.rm = TRUE)
    data.frame(
      treatment = trt,
      x     = as.numeric(factor(trt, levels = treatment_levels)),
      y     =3 ,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  })
  bind_rows(result)
}
 

###### Exons

###### Define tables to plot
treatments_keep <- "Speckle-dTAGvsNEG_6h_exonsOnly"
#"Speckle-dTAGvsNEG_6h_exonsOnly"
combined_sub <- DGE_exonTCReads_combined %>% filter(treatment %in% treatments_keep)
sig_only_sub <- combined_sub %>% filter(Significance == "Significant")

treatment_order <- c("Speckle-dTAGvsNEG_6h_exonsOnly")#"Speckle-dTAGvsNEG_6h_exonsOnly")#  # set your desired order here

combined_sub$treatment <- factor(combined_sub$treatment, levels = treatment_order)
sig_only_sub$treatment <- factor(sig_only_sub$treatment, levels = treatment_order)

treatment_levels <- levels(combined_sub$treatment)
treatment_levels <- levels(factor(combined_sub$treatment))

## Helper functions 
compute_annotations <- function(all_df, group_col = "treatment") {
  treatments <- unique(all_df[[group_col]])
  result <- lapply(treatments, function(trt) {
    all_vals <- all_df$Log2FC[all_df[[group_col]] == trt & !is.na(all_df$Log2FC)]
    n_total  <- length(all_vals)
    n_up     <- sum(all_df[[group_col]] == trt & all_df$Color_factor == "Significantly_Up",   na.rm = TRUE)
    n_down   <- sum(all_df[[group_col]] == trt & all_df$Color_factor == "Significantly_Down", na.rm = TRUE)
    data.frame(
      treatment = trt,
      x     = as.numeric(factor(trt, levels = treatment_levels)),
      y     =3 ,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  })
  bind_rows(result)
}

### Compute plot data
scale_factor <- 0.2

left_dens  <- compute_half_violin(combined_sub, "left")
right_dens <- compute_half_violin(sig_only_sub, "right")

left_dens$x_pos  <- as.numeric(factor(left_dens$treatment,  levels = treatment_levels)) + left_dens$density  * scale_factor
right_dens$x_pos <- as.numeric(factor(right_dens$treatment, levels = treatment_levels)) + right_dens$density * scale_factor

box_left  <- compute_boxplot_stats(combined_sub, "left")
box_right <- compute_boxplot_stats(sig_only_sub, "right")
box_width <- 0.04

annot_df <- compute_annotations(combined_sub)

### Plot 
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")
pdf("SpdTAGvsNeg_SplitViolin_SLAMseq_exons.pdf", height = 5, width = 7)
ggplot() +
  geom_polygon(data = left_dens,  aes(x = x_pos, y = y, group = treatment, fill = "All"),         alpha = 0.7, color = "black", linewidth = 0.8) +
  geom_polygon(data = right_dens, aes(x = x_pos, y = y, group = treatment, fill = "Significant"), alpha = 0.7, color = "black", linewidth = 0.8) +
  geom_boxplot(data = box_left,  aes(x = x, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, group = treatment),
               stat = "identity", width = box_width, fill = "white", alpha = 0.8) +
  geom_boxplot(data = box_right, aes(x = x, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, group = treatment),
               stat = "identity", width = box_width, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = annot_df, aes(x = x, y = y, label = label),
            vjust = -0.3, size = 8, lineheight = 0.9) +
  scale_x_continuous(breaks = seq_along(treatment_levels), 
                     labels = treatment_levels,
                     limits = c(0.7, 1.3)) +
  scale_fill_manual(values = c("All" = "#7fbfff", "Significant" = "#ff7f7f")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme_classic(base_size = 14) +   
  theme(
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.ticks.x  = element_line(linewidth = 2.7, lineend = "round"), 
    axis.ticks.y  = element_line(linewidth = 2.7, lineend = "round"), 
    axis.ticks.length = unit(0.3, "cm"), 
    axis.line     = element_line(linewidth = 2.7, lineend = "round"),  
    axis.text.y = element_text(size = 37, margin = margin(r = 5)),
    axis.title.y  = element_blank(), 
    legend.position = "none"
  ) +
  labs(y = "Log2 Fold Change")
dev.off()

##########


############## Make barplots for intron overrepresentation by SPIN for U1 and Speckle-dTAG  ####################


##############  Subset genes with intron upregulation in SP-dTAG treatments

gene_ids <- DGE_intronTCReads_combined[DGE_intronTCReads_combined$treatment == "Speckle-dTAGvsNEG_6h_intronssOnly" & DGE_intronTCReads_combined$Color_factor == "Significantly_Up",]

gene_ids$Gene <- sub("_[^_]*$", "", gene_ids$Gene)

gene_ids <- as.vector(gene_ids$Gene)


########  Fetch gene coordinates from Ensembl 

fetch_ensembl_coords <- function(ids, batch_size = 1000) {
  base_url <- "https://rest.ensembl.org/lookup/id"
  results  <- list()
  
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  
  for (i in seq_along(batches)) {
    message("  Querying Ensembl batch ", i, "/", length(batches), " ...")
    body <- toJSON(list(ids = batches[[i]]), auto_unbox = FALSE)
    
    resp <- POST(
      url    = base_url,
      body   = body,
      encode = "raw",
      add_headers(
        "Content-Type" = "application/json",
        "Accept"       = "application/json"
      )
    )
    
    if (http_error(resp)) {
      warning("Ensembl API error on batch ", i, ": ", status_code(resp))
      next
    }
    
    parsed <- fromJSON(rawToChar(resp$content), simplifyVector = FALSE)
    results <- c(results, parsed)
    Sys.sleep(0.5)  
  }
  
  # Flatten into a data frame
  coords <- lapply(names(results), function(id) {
    g <- results[[id]]
    if (is.null(g) || is.null(g$seq_region_name)) return(NULL)
    data.frame(
      gene_id    = id,
      gene_name  = if (!is.null(g$display_name)) g$display_name else NA_character_,
      chrom      = paste0("chr", g$seq_region_name),
      start      = as.integer(g$start),
      end        = as.integer(g$end),
      strand     = g$strand,
      biotype    = if (!is.null(g$biotype)) g$biotype else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(Filter(Negate(is.null), coords))
}

gene_coords <- fetch_ensembl_coords(gene_ids)

# Filter chromosomes
canonical <- paste0("chr", c(1:22, "X"))
gene_coords_filtered <- gene_coords %>%
  mutate(TSS = ifelse(strand == 1, start, end))
# Save gene coordinate 
gene_bed <- gene_coords_filtered %>%
  transmute(chrom, start = start - 1L, end, name = gene_id, score = 0, strand)
bin_50kb <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Sequencing/data/binning/hg38_50kb.bed",
                       sep="\t", header=FALSE)
colnames(bin_50kb) <- c("chrom","start","end")

bin_150kb <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Sequencing/data/binning/hg38_150kb.bed",
                       sep="\t", header=FALSE)
colnames(bin_150kb) <- c("chrom","start","end")
# Filter chromosomes
bin_50kb <- bin_50kb %>%
  filter(chrom %in% canonical) %>%
  arrange(chrom, start)

bin_150kb <- bin_150kb %>%
  filter(chrom %in% canonical) %>%
  arrange(chrom, start)

gene_gr <- with(gene_coords_filtered,
                GRanges(seqnames = chrom,
                        ranges   = IRanges(start = TSS, end = TSS))) 

bin_gr <- with(bin_50kb,
               GRanges(seqnames = chrom,
                       ranges   = IRanges(start = start + 1L, end = end)))

bin_gr150 <- with(bin_150kb,
               GRanges(seqnames = chrom,
                       ranges   = IRanges(start = start + 1L, end = end)))

hits <- findOverlaps(gene_gr, bin_gr)
hits150 <- findOverlaps(gene_gr, bin_gr150)

gene_counts <- data.frame(
  bin_idx  = subjectHits(hits),
  gene_idx = queryHits(hits)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_SpdTAG = n_distinct(gene_idx), .groups = "drop")

gene_counts150 <- data.frame(
  bin_idx  = subjectHits(hits150),
  gene_idx = queryHits(hits150)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_SpdTAG = n_distinct(gene_idx), .groups = "drop")

# Attach counts to bins (0 for empty bins)
bin_50kb$bin_idx   <- seq_len(nrow(bin_50kb))
coverage_bed <- bin_50kb %>%
  left_join(gene_counts, by = "bin_idx") %>%
  mutate(gene_count_SpdTAG = replace(gene_count_SpdTAG, is.na(gene_count_SpdTAG), 0L)) %>%
  select(chrom, start, end, gene_count_SpdTAG)

bin_150kb$bin_idx   <- seq_len(nrow(bin_150kb))
coverage_bed150 <- bin_150kb %>%
  left_join(gene_counts150, by = "bin_idx") %>%
  mutate(gene_count_SpdTAG = replace(gene_count_SpdTAG, is.na(gene_count_SpdTAG), 0L)) %>%
  select(chrom, start, end, gene_count_SpdTAG)



# Merge on chrom + start + end (exact bin match)
df_residuals <- merge(df_residuals, coverage_bed, by=c("chrom", "start", "end"), all.x= TRUE)
df_residuals_150 <- merge(df_residuals_150, coverage_bed150, by=c("chrom", "start", "end"), all.x= TRUE)



######
##############  Subset genes downregulated in SP-dTAG treatments

gene_ids <- DGE_allTCReads_combined[DGE_allTCReads_combined$treatment == "Speckle-dTAGvsNEG_6h_allReads" & DGE_allTCReads_combined$Color_factor == "Significantly_Down",]

gene_ids$Gene <- sub("_[^_]*$", "", gene_ids$Gene)

gene_ids <- as.vector(gene_ids$Gene)

########  Fetch gene coordinates from Ensembl 
gene_coords <- fetch_ensembl_coords(gene_ids)

# Filter chromosomes
canonical <- paste0("chr", c(1:22, "X"))
gene_coords_filtered <- gene_coords %>%
  mutate(TSS = ifelse(strand == 1, start, end))
# Save gene coordinate 
gene_bed <- gene_coords_filtered %>%
  transmute(chrom, start = start - 1L, end, name = gene_id, score = 0, strand)

# Filter chromosomes

gene_gr <- with(gene_coords_filtered,
                GRanges(seqnames = chrom,
                        ranges   = IRanges(start = TSS, end = TSS))) 

hits <- findOverlaps(gene_gr, bin_gr)
hits150 <- findOverlaps(gene_gr, bin_gr150)

gene_counts <- data.frame(
  bin_idx  = subjectHits(hits),
  gene_idx = queryHits(hits)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_SpdTAG_down = n_distinct(gene_idx), .groups = "drop")

gene_counts150 <- data.frame(
  bin_idx  = subjectHits(hits150),
  gene_idx = queryHits(hits150)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_SpdTAG_down = n_distinct(gene_idx), .groups = "drop")

# Attach counts to bins (0 for empty bins)

coverage_bed <- bin_50kb %>%
  left_join(gene_counts, by = "bin_idx") %>%
  mutate(gene_count_SpdTAG_down  = replace(gene_count_SpdTAG_down , is.na(gene_count_SpdTAG_down), 0L)) %>%
  select(chrom, start, end, gene_count_SpdTAG_down)

coverage_bed150 <- bin_150kb %>%
  left_join(gene_counts150, by = "bin_idx") %>%
  mutate(gene_count_SpdTAG_down  = replace(gene_count_SpdTAG_down , is.na(gene_count_SpdTAG_down), 0L)) %>%
  select(chrom, start, end, gene_count_SpdTAG_down)



# Merge on chrom + start + end (exact bin match)
df_residuals <- merge(df_residuals, coverage_bed, by=c("chrom", "start", "end"), all.x= TRUE)
df_residuals_150 <- merge(df_residuals_150, coverage_bed150, by=c("chrom", "start", "end"), all.x= TRUE)






#####
##############  Subset genes with intron upregulation in U1 AMO treatments

gene_ids <- DGE_intronTCReads_combined[DGE_intronTCReads_combined$treatment == "U1AMO_4h_30uM_intronsOnly" & DGE_intronTCReads_combined$Color_factor == "Significantly_Up",]

gene_ids$Gene <- sub("_[^_]*$", "", gene_ids$Gene)

gene_ids <- as.vector(gene_ids$Gene)


########  Fetch gene coordinates from Ensembl 

gene_coords <- fetch_ensembl_coords(gene_ids)

# Filter chromosomes
canonical <- paste0("chr", c(1:22, "X"))
gene_coords_filtered <- gene_coords %>%
  mutate(TSS = ifelse(strand == 1, start, end))
# Save gene coordinate 
gene_bed <- gene_coords_filtered %>%
  transmute(chrom, start = start - 1L, end, name = gene_id, score = 0, strand)

# Filter chromosomes
gene_gr <- with(gene_coords_filtered,
                GRanges(seqnames = chrom,
                        ranges   = IRanges(start = TSS, end = TSS)))  

hits <- findOverlaps(gene_gr, bin_gr)
hits150 <- findOverlaps(gene_gr, bin_gr150)

gene_counts <- data.frame(
  bin_idx  = subjectHits(hits),
  gene_idx = queryHits(hits)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_U1 = n_distinct(gene_idx), .groups = "drop")

gene_counts150 <- data.frame(
  bin_idx  = subjectHits(hits150),
  gene_idx = queryHits(hits150)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_U1 = n_distinct(gene_idx), .groups = "drop")

# Attach counts to bins (0 for empty bins)
coverage_bed <- bin_50kb %>%
  left_join(gene_counts, by = "bin_idx") %>%
  mutate(gene_count_U1 = replace(gene_count_U1, is.na(gene_count_U1), 0L)) %>%
  select(chrom, start, end, gene_count_U1)

coverage_bed150 <- bin_150kb %>%
  left_join(gene_counts150, by = "bin_idx") %>%
  mutate(gene_count_U1 = replace(gene_count_U1, is.na(gene_count_U1), 0L)) %>%
  select(chrom, start, end, gene_count_U1)


# Merge on chrom + start + end (exact bin match)
df_residuals <- merge(df_residuals, coverage_bed, by=c("chrom", "start", "end"), all.x= TRUE)
df_residuals_150 <- merge(df_residuals_150, coverage_bed150, by=c("chrom", "start", "end"), all.x= TRUE)

######
##############  Subset genes downregulated in U14h treatments

gene_ids <- DGE_allTCReads_combined[DGE_allTCReads_combined$treatment == "U1AMO_4h_30uM_allReads" & DGE_allTCReads_combined$Color_factor == "Significantly_Down",]

gene_ids$Gene <- sub("_[^_]*$", "", gene_ids$Gene)

gene_ids <- as.vector(gene_ids$Gene)

########  Fetch gene coordinates from Ensembl 
gene_coords <- fetch_ensembl_coords(gene_ids)

# Filter chromosomes
canonical <- paste0("chr", c(1:22, "X"))
gene_coords_filtered <- gene_coords %>%
  mutate(TSS = ifelse(strand == 1, start, end))
# Save gene coordinate 
gene_bed <- gene_coords_filtered %>%
  transmute(chrom, start = start - 1L, end, name = gene_id, score = 0, strand)

# Filter chromosomes

gene_gr <- with(gene_coords_filtered,
                GRanges(seqnames = chrom,
                        ranges   = IRanges(start = TSS, end = TSS)))  

hits <- findOverlaps(gene_gr, bin_gr)
hits150 <- findOverlaps(gene_gr, bin_gr150)

gene_counts <- data.frame(
  bin_idx  = subjectHits(hits),
  gene_idx = queryHits(hits)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_U14h_down = n_distinct(gene_idx), .groups = "drop")

gene_counts150 <- data.frame(
  bin_idx  = subjectHits(hits150),
  gene_idx = queryHits(hits150)
) %>%
  group_by(bin_idx) %>%
  summarise(gene_count_U14h_down = n_distinct(gene_idx), .groups = "drop")

# Attach counts to bins (0 for empty bins)
coverage_bed <- bin_50kb %>%
  left_join(gene_counts, by = "bin_idx") %>%
  mutate(gene_count_U14h_down  = replace(gene_count_U14h_down , is.na(gene_count_U14h_down), 0L)) %>%
  select(chrom, start, end, gene_count_U14h_down)

coverage_bed150 <- bin_150kb %>%
  left_join(gene_counts150, by = "bin_idx") %>%
  mutate(gene_count_U14h_down  = replace(gene_count_U14h_down , is.na(gene_count_U14h_down), 0L)) %>%
  select(chrom, start, end, gene_count_U14h_down)


# Merge on chrom + start + end (exact bin match)
df_residuals <- merge(df_residuals, coverage_bed, by=c("chrom", "start", "end"), all.x= TRUE)
df_residuals_150 <- merge(df_residuals_150, coverage_bed150, by=c("chrom", "start", "end"), all.x= TRUE)

df_residuals$chrom <- factor(df_residuals$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals <- df_residuals[order(df_residuals$chrom, df_residuals$start,df_residuals$end),]

df_residuals_150$chrom <- factor(df_residuals_150$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals_150 <- df_residuals_150[order(df_residuals_150$chrom, df_residuals_150$start,df_residuals_150$end),]


b <- df_residuals_150
b <- b[,c("chrom","start", "end","SPIN","seg_1","Speckle_seg","Interior_Act1_seg","Interior_Act2_seg","Interior_Repr1_seg", "Interior_Repr2_seg","Near_Lm1_seg","Near_Lm2_seg",
          "Lamina_seg","gene_count_SpdTAG","gene_count_SpdTAG_down","gene_count_U1","gene_count_U14h_down")]


b$SpdTAG_intUP_seg <- 1
b$SpdTAG_intUP_seg <- ifelse(b$gene_count_SpdTAG > 0, b$SpdTAG_intUP_seg, 10)

b$chrom <- factor(b$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
b <- b[order(b$chrom, b$start,b$end),]



#######  Fig. 4f

##### Make segment cols U1

b$U1_intUP_seg <- 3
b$U1_intUP_seg <- ifelse(b$gene_count_U1 > 0, b$U1_intUP_seg, 10)


b$chrom <- factor(b$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
b <- b[order(b$chrom, b$start,b$end),]


seg_size <- 24
lw       <- 6
ft       <- 2.6

chroms_to_plot <- c("chr1", "chr3", "chr6", "chr11", "chr16", "chr18", "chr19")

# global x-limit driven by the longest chromosome in the selection
global_xlim <- c(0, max(b$end[b$chrom %in% chroms_to_plot]) / 1e6)
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")
pdf("SPINsegs_IntronHitsegs_multichr7_150kb_alt.pdf", height = length(chroms_to_plot) * 1.2, width = 10)

heights <- c(rep(1, length(chroms_to_plot) - 1), 1.3) 
layout(matrix(1:length(chroms_to_plot), ncol = 1), heights = heights)

for (i in chroms_to_plot) {
  is_last <- i == chroms_to_plot[length(chroms_to_plot)]
  par(mar = c(if (is_last) 1.5 else 0.1, 4, 0.1, 2))
  
  x <- b[b$chrom == i, ] 
  
  plot(x$start / 1e6, x$seg_1, type = "n",
       xlim = global_xlim, xlab = "", ylab = "",
       axes = FALSE, ylim = c(-3.5, 4))
  
  # SPIN state
  segments(x0=x$end/1e6, y0=x$Speckle_seg,        x1=x$start/1e6, y1=x$Speckle_seg,        lend=1, lwd=seg_size, col="#9E0142")
  segments(x0=x$end/1e6, y0=x$Interior_Act1_seg,  x1=x$start/1e6, y1=x$Interior_Act1_seg,  lend=1, lwd=seg_size, col="#F46D43")
  segments(x0=x$end/1e6, y0=x$Interior_Act2_seg,  x1=x$start/1e6, y1=x$Interior_Act2_seg,  lend=1, lwd=seg_size, col="#FDAE61")
  segments(x0=x$end/1e6, y0=x$Interior_Repr1_seg, x1=x$start/1e6, y1=x$Interior_Repr1_seg, lend=1, lwd=seg_size, col="#FFFFBF")
  segments(x0=x$end/1e6, y0=x$Interior_Repr2_seg, x1=x$start/1e6, y1=x$Interior_Repr2_seg, lend=1, lwd=seg_size, col="#ABDDA4")
  segments(x0=x$end/1e6, y0=x$Near_Lm1_seg,       x1=x$start/1e6, y1=x$Near_Lm1_seg,       lend=1, lwd=seg_size, col="#66C2A5")
  segments(x0=x$end/1e6, y0=x$Near_Lm2_seg,       x1=x$start/1e6, y1=x$Near_Lm2_seg,       lend=1, lwd=seg_size, col="#3288BD")
  segments(x0=x$end/1e6, y0=x$Lamina_seg,          x1=x$start/1e6, y1=x$Lamina_seg,          lend=1, lwd=seg_size, col="#5E4FA2")
  
  # SpdTAG intron UP
  segments(x0=x$end/1e6, y0=x$SpdTAG_intUP_seg, x1=x$start/1e6, y1=x$SpdTAG_intUP_seg,
           lend=1, lwd=seg_size, col="black")
  
  # U1 intron UP
  segments(x0=x$end/1e6, y0=x$U1_intUP_seg, x1=x$start/1e6, y1=x$U1_intUP_seg,
           lend=1, lwd=seg_size, col="black")
  
  if (is_last) axis(1, lwd = lw, labels = FALSE)
}

dev.off()






####### Fig. 4d


lw <- 0.9

NT <- df_residuals[df_residuals$SPIN %in% c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),]

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")
###### Count genes with intron retention by SPIN state ####
pdf("Intron_overRepresentation_bySPIN_SpdTAG_alt.pdf", height = 1, width = 4)
NT %>%
  filter(!is.na(SPIN), !is.na(gene_count_SpdTAG)) %>%
  group_by(SPIN) %>%
  summarise(total = sum(gene_count_SpdTAG)) %>%
  ggplot(aes(x = SPIN, y = total, fill = SPIN)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  theme(axis.ticks = element_line(colour = "black", linewidth = lw, lineend = "round")) +
  theme(axis.ticks.length = unit(0.2, "cm")) +
  theme(axis.line = element_line(linewidth = lw, color = "black", lineend = "round")) +
  labs(x = NULL, y = NULL)
dev.off()



pdf("Intron_overRepresentation_bySPIN_U1AMO_alt.pdf", height = 1, width = 4)
NT %>%
  filter(!is.na(SPIN), !is.na(gene_count_U1)) %>%
  group_by(SPIN) %>%
  summarise(total = sum(gene_count_U1)) %>%
  ggplot(aes(x = SPIN, y = total, fill = SPIN)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  theme(axis.ticks = element_line(colour = "black", linewidth = lw, lineend = "round")) +
  theme(axis.ticks.length = unit(0.2, "cm")) +
  theme(axis.line = element_line(linewidth = lw, color = "black", lineend = "round")) +
  labs(x = NULL, y = NULL)
dev.off()



####### Fig. 4e

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")
###### Count genes with intron retention by SPIN state ####
pdf("DownregulatedGene_bySPIN_SpdTAG_alt.pdf", height = 1, width = 4)
NT %>%
  filter(!is.na(SPIN), !is.na(gene_count_SpdTAG_down)) %>%
  group_by(SPIN) %>%
  summarise(total = sum(gene_count_SpdTAG_down)) %>%
  ggplot(aes(x = SPIN, y = total, fill = SPIN)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  theme(axis.ticks = element_line(colour = "black", linewidth = lw, lineend = "round")) +
  theme(axis.ticks.length = unit(0.2, "cm")) +
  theme(axis.line = element_line(linewidth = lw, color = "black", lineend = "round")) +
  labs(x = NULL, y = NULL)
dev.off()


pdf("DownregulatedGenes_bySPIN_U1AMO_alt.pdf", height = 1, width = 4)
NT %>%
  filter(!is.na(SPIN), !is.na(gene_count_U14h_down)) %>%
  group_by(SPIN) %>%
  summarise(total = sum(gene_count_U14h_down)) %>%
  ggplot(aes(x = SPIN, y = total, fill = SPIN)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_blank(),
    legend.position = "none"
  ) +
  theme(axis.ticks = element_line(colour = "black", linewidth = lw, lineend = "round")) +
  theme(axis.ticks.length = unit(0.2, "cm")) +
  theme(axis.line = element_line(linewidth = lw, color = "black", lineend = "round")) +
  labs(x = NULL, y = NULL)
dev.off()




##### Fig. 4g-k

###################   Tracks  #########################

seg_size <- 30

begin <- 0
term <- 80000000

t <- 1.6

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/tracks")

x <- df_residuals[df_residuals$chrom == "chr11",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_TREATdif_track_ch11_dTAG_wideA_SPINonly_alt.pdf", height = 1.6, width = 3.7)
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


 axis(1, lwd=0.4, cex.axis=t, tcl=-0.3, labels=FALSE) 

dev.off()




begin <- 0
term <- 80000000
t <- 2
lw <- 2.3
  
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/tracks")

x <- df_residuals_150[df_residuals_150$chrom == "chr11",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_TREATdif_track_ch11_dTAG_wideB_alt.pdf", height = 7, width = 3.7)
par(mfrow=c(4,1), mar=c(1.5, 4, 2, 2) + 0.1,oma=c(5, 0, 0, 0), mgp=c(3, 1, 0))
plot(x$start/1000000, x$LOS_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif, type="l", col="black", 
     ylim=c(-0.3, 0.22), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw, cex.axis=t, labels=FALSE) 
axis(2, lwd=lw, cex.axis=t)
box(bty="l", lwd=lw)
  
plot(x$start/1000000, x$LOS_residuals_range2Mb_K562SpeckledTAG_11hdTAG120mPD_DpnII_R1_20240828_dif, type="l", col="black",
     ylim=c(-0.3, 0.22), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw, cex.axis=t, labels=FALSE)
axis(2, lwd=lw, cex.axis=t)
box(bty="l", lwd=lw)
  
plot(x$start/1000000, x$LOS_residuals_range2Mb_K562SpeckledTAG_NegNTTD120mPD_DpnII_R1_20240828_dif, type="l", col="black",
     ylim=c(-0.3, 0.22), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw, cex.axis=t, labels=FALSE)
axis(2, lwd=lw, cex.axis=t)
box(bty="l", lwd=lw)
  
  
plot(x$start/1000000, x$LOS_residuals_range2Mb_K562SpeckledTAG_dTAGNTTD120mPD_DpnII_R1_20240828_dif, type="l", col="black", 
     ylim=c(-0.3, 0.22), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw, cex.axis=t, labels=TRUE, mgp=c(3, 1.2, 0)) 
axis(2, lwd=lw, cex.axis=t)
box(bty="l", lwd=lw)
  
dev.off()



################  Global violins    ########################
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/violin")

NT <- df_residuals_150[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif','LOS_residuals_range2Mb_K562SpeckledTAG_11hdTAG120mPD_DpnII_R1_20240828_dif','LOS_residuals_range2Mb_K562SpeckledTAG_NegNTTD120mPD_DpnII_R1_20240828_dif','LOS_residuals_range2Mb_K562SpeckledTAG_dTAGNTTD120mPD_DpnII_R1_20240828_dif')]

NTTD <- NT[,c(4, 6)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "dTAG_6h"
NTTD <- NTTD[NTTD$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

lw <- 1.4

#Violin plotting the LOS residuals for dTAG 6h
pdf("LOSr_SPIN_Speckle_IntRepr2_by_dTAG6h_violin_alt.pdf", width=4.5, height=2.05)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#ABDDA4",  "#5E4FA2")) +
  ylim(-0.295, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = lw, lineend = "round"),  
        axis.ticks = element_line(linewidth = lw, lineend = "round"), 
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

options(scipen = 999)
aggregate(LOS_residuals ~ SPIN, data = na.omit(NTTD), FUN = mean)


U14h <- NT[,c(4, 7)]
colnames(U14h)[2] <- "LOS_residuals"
U14h$treat <- "dTAG_11h"
U14h <- U14h[U14h$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/violin")
#Violin plotting the LOS residuals for dTag10h
pdf("LOSr_SPIN_Speckle_IntRepr2_by_dTag11h_violin_alt.pdf", width=4.5, height=2.05)
ggplot(na.omit(U14h), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#ABDDA4",  "#5E4FA2")) +
  ylim(-0.295, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2) +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = lw, lineend = "round"),  
        axis.ticks = element_line(linewidth = lw, lineend = "round"), 
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

aggregate(LOS_residuals ~ SPIN, data = na.omit(U14h), FUN = mean)


U18h <- NT[,c(4, 8)]
colnames(U18h)[2] <- "LOS_residuals"
U18h$treat <- "TPL/DRB"
U18h <- U18h[U18h$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/violin")
#Violin plotting the LOS residuals for TD
pdf("LOSr_SPIN_Speckle_IntRepr2_by_TD_violin_alt.pdf", width=4.5, height=2.05)
ggplot(na.omit(U18h), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#ABDDA4",  "#5E4FA2")) +
  ylim(-0.295, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2) +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = lw, lineend = "round"), 
        axis.ticks = element_line(linewidth = lw, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

aggregate(LOS_residuals ~ SPIN, data = na.omit(U18h), FUN = mean)


RIRA <- NT[,c(4, 9)]
colnames(RIRA)[2] <- "LOS_residuals"
RIRA$treat <- "dTAG/TPL/DRB"
RIRA <- RIRA[RIRA$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/violin")
#Violin plotting the LOS residuals for dTag6h + TD
pdf("LOSr_SPIN_Speckle_IntRepr2_by_dTag6hTD_violin_alt.pdf", width=4.5, height=2.05)
ggplot(na.omit(RIRA), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#ABDDA4",  "#5E4FA2")) +
  ylim(-0.295, 0.22) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2) +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = lw, lineend = "round"),  
        axis.ticks = element_line(linewidth = lw, lineend = "round"),   
        axis.ticks.length = unit(0.29, "cm"))

dev.off()

aggregate(LOS_residuals ~ SPIN, data = na.omit(RIRA), FUN = mean)



####### Extended Fig. 4c
######### Crane maps were generated using "make_cranemaps_EDFig4.ipynb"

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")
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

#####



## Fig.4c - LOS residuals Speckle degron cell line
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")

x <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/150kb/LOS_residuals/LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828.bedGraph",
                sep="\t", header=TRUE)
x <-x[x$chrom == "chr11", ]
pdf("LOSres_NegDMSO_track_ch11_xaxis_150kb.pdf", height = 2.6, width = 9)
par(mfrow=c(1,1), mar=c(3, 4, 2, 2) + 0.1,mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=TRUE) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()



######## Extended data Fig. 4d -  Global violins    ########################
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_speckle_degron/")

NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828')]


NTTD <- NT[,c(4, 6)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "dTAG10h"


NTTD <- NTTD[NTTD$SPIN == c("Speckle", "Interior_Act1", "Interior_Act2","Interior_Repr1", "Interior_Repr2", "Near_Lm1", "Near_Lm2", "Lamina"),]

#Violin plotting the LOS residuals for dTAG 6h
pdf("LOSr_allSPIN_global_dTAGDMSO_violin_alt.pdf", width=10, height=3.5)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.32, 0.22) +
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


