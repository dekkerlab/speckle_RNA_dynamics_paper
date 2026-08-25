library(tidyverse)
library(ggbeeswarm)



#### df_residuals table was created using the script "C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/scripts/make_residual_table_toggleRES_hg38_20260225.R"
#### cis_percent .bed files were created using ""

df_residuals <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/dataframes/df_residuals_20260604_50kb.bed",
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



setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/violin")


df_residuals_150kb <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/150kb/dataframes/df_residuals_20260604_150kb.bed", header=TRUE, sep="\t")
df_residuals_150kb$SPIN <- factor(df_residuals_150kb$SPIN, order = TRUE, levels =c("Speckle", "Interior_Act1", "Interior_Act2", "Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina","Lamina_Like"))

df_residuals_150kb <- df_residuals_150kb %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),
                  ordered = TRUE))

################################ imaging analsysis Figs3b and c #########################################################

##### Aivia output was munged using "C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/scripts/munge_AiviaOutput.R" and the master_table was generated)

master_table <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Microscopy Analysis/20260129d_FigureRNAtreat/2026-01-29-17-09-34 IF_analysis/combined_measurements_2026-01-29-16-30-54 IF_analysis_d.csv",
                           sep=",", header=TRUE)


##########  Add treatment columns 
master_table$treatment <- NA

master_table$treatment[grepl("DMSO", master_table$Image)] <- "DMSO"
master_table$treatment[grepl("TD", master_table$Image)] <- "TD"
master_table$treatment[grepl("cAMO", master_table$Image)] <- "cAMO"
master_table$treatment[grepl("contAMO", master_table$Image)] <- "cAMO"
master_table$treatment[grepl("U1AMO", master_table$Image)] <- "U1AMO"


###### number images from each treatment

image_number <- NA

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("DMSO", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("DMSO", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("TD", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("TD", treatment)]))),
    image_number  
  ))


master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("cAMO", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("cAMO", treatment)]))),
    image_number  
  ))

master_table <- master_table %>%
  mutate(image_number = if_else(
    grepl("U1AMO", treatment),
    as.integer(factor(Image, levels = unique(Image[grepl("U1AMO", treatment)]))),
    image_number  
  ))



###### Sort table

master_table <- master_table %>%
  mutate(treatment = factor(treatment, 
                            levels = c("DMSO", "TD", "cAMO","U1AMO")))


master_table <- master_table %>%
  arrange(treatment, image_number)


###### Add identifier columns

master_table$sample <- NA
master_table$sample[master_table$treatment== "DMSO"] <- "control"
master_table$sample[master_table$treatment== "TD"] <- "treatment"
master_table$sample[master_table$treatment== "cAMO"] <- "control"
master_table$sample[master_table$treatment== "U1AMO"] <- "treatment"


master_table$treat_id <- NA
master_table$treat_id[master_table$treatment== "DMSO"] <- "TPL/DRB"
master_table$treat_id[master_table$treatment== "TD"] <- "TPL/DRB"
master_table$treat_id[master_table$treatment== "cAMO"] <- "U1 AMO"
master_table$treat_id[master_table$treatment== "U1AMO"] <- "U1 AMO"


#### Convert pixel values to um

#### Voxel Size
x <- 0.038449
y <- 0.038449
z <- 0.294583333

master_table <- master_table %>%
  mutate(Volume..µm.. = Volume..px.. * x * y * z) %>%
  relocate(Volume..µm.., .after = 3)

long_vol <- master_table[,c('Image','treat_id','sample','Volume..µm..')]
long_vol <- long_vol[long_vol$treat_id ==  c("U1 AMO","TPL/DRB"),]


# Calculate counts for each dodged group
count_data <- long_vol %>%
  group_by(treat_id, sample) %>%
  summarise(
    n = n(),
    y_pos = 35,
    .groups = 'drop'
  )
long_vol$group_color <- paste(long_vol$treat_id, long_vol$sample, sep = ".")


## Fig.2b - Speckle volume by treatment (IF)

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
pdf("Speckle_vol_by_treatment_IF_v10.pdf", height = 5, width = 6.5)
ggplot(na.omit(long_vol), aes(treat_id, Volume..µm.., color = group_color, group = interaction(sample, treat_id))) +
  geom_violin(aes(group = interaction(sample, treat_id), fill = group_color), 
              position = position_dodge(width = 0.75), color = NA) +
  geom_boxplot(aes(alpha = NULL, color = NULL), outlier.colour = NA, width = 0.15, 
               position = position_dodge(width = 0.75), fill = "white") +
  geom_text(data = count_data, aes(x = treat_id, y = y_pos, label = n, group = interaction(sample, treat_id), color = NULL),
            position = position_dodge(width = 0.75), size = 9, vjust = -0.5, 
            show.legend = FALSE, color = "black", alpha = 1) +
  scale_y_log10(limits = c(0.05, 60)) +
  scale_fill_manual(values = c(
    "TPL/DRB.control" = "red",
    "TPL/DRB.treatment" = "blue",
    "U1 AMO.control" = "red",
    "U1 AMO.treatment" = "blue"
  )) +
  scale_color_manual(values = c(
    "TPL/DRB.control" = "red",
    "TPL/DRB.treatment" = "blue",
    "U1 AMO.control" = "red",
    "U1 AMO.treatment" = "blue"
  )) +
  theme_minimal(base_size = 20) +
  labs(x = NULL) +
  ylab(bquote('Volume '(µm^3))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 2.3, lineend = "round"),
        axis.ticks.length = unit(0.3, "cm"),
        axis.line = element_line(colour = "black", linewidth = 2.3, lineend = "round"),
        axis.text = element_text(size = 40, colour = "black"), 
        axis.title = element_text(size = 30, face = "bold"),
        legend.text = element_text(size = 25), 
        legend.title = element_text(size = 30),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

dev.off()


# Filter for control samples only
control_data <- long_vol %>%
  filter(sample == "control")

# Get the two treatment groups
treat1_data <- control_data$Volume..µm..[control_data$treat_id == "TPL/DRB"]
treat2_data <- control_data$Volume..µm..[control_data$treat_id == "U1 AMO"]  

# Perform KS test
ks_result <- ks.test(treat1_data, treat2_data)

# Print results
cat("KS test between control samples:\n")
cat("TPL/DRB vs", unique(control_data$treat_id)[2], "\n")
cat("D statistic:", ks_result$statistic, "\n")
cat("p-value:", ks_result$p.value, "\n")


# Perform KS tests for each treatment
ks_results <- long_vol %>%
  group_by(treat_id) %>%
  summarise(
    ks_statistic = ks.test(Volume..µm..[sample == "control"], 
                           Volume..µm..[sample == "treatment"])$statistic,
    ks_pvalue = ks.test(Volume..µm..[sample == "control"], 
                        Volume..µm..[sample == "treatment"])$p.value,
    .groups = "drop"
  )

print(ks_results)



######## Plot number of speckles per cell
w <- as.data.frame(table(master_table$Image))

for (i in c(1:length(w$Var1))) {
  w$treat_id[i] <- paste(unlist(str_split(w$Var1[i],"_"))[4])
}
w$treat_id <- gsub("contAMO", "cAMO", w$treat_id)
w$treat_id <- factor(w$treat_id, order = TRUE, levels =c("DMSO","TD", "cAMO","U1AMO")) 

w$sample <- NA
w$sample[w$treat_id== "DMSO"] <- "control"
w$sample[w$treat_id== "TD"] <- "treatment"
w$sample[w$treat_id== "cAMO"] <- "control"
w$sample[w$treat_id== "U1AMO"] <- "treatment"

w$treatment <- NA
w$treatment[w$treat_id== "DMSO"] <- "TPL/DRB"
w$treatment[w$treat_id== "TD"] <- "TPL/DRB"
w$treatment[w$treat_id== "cAMO"] <- "U1 AMO"
w$treatment[w$treat_id== "U1AMO"] <- "U1 AMO"


# Calculate counts for each dodged group
count_data <- w %>%
  group_by(treatment, sample) %>%
  summarise(
    n = n(),
    y_pos = 106,  
    .groups = 'drop'
  )

w <- w %>%
  mutate(group_color = interaction(treatment, sample))

## Fig.2c - Speckle count by treatment (IF)

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
pdf("Speckle_count_by_treatment_IF_v10.pdf", height = 5, width = 6.1)
ggplot(na.omit(w), aes(x = treatment, y = Freq, color = group_color, 
                       group = interaction(sample, treatment))) +
  geom_boxplot(aes(alpha = NULL, color = NULL), outlier.colour = NA, width = 0.5,
               position = position_dodge(width = 0.75), fill = "white") +
  geom_beeswarm(aes(group = interaction(sample, treatment)), 
                size = 3, dodge.width = 0.75) +
  geom_text(data = count_data, 
            aes(x = treatment, y = y_pos, label = n, 
                group = interaction(sample, treatment), color = NULL),
            position = position_dodge(width = 0.75),
            size = 9, vjust = -0.5, show.legend = FALSE,
            color = "black", alpha = 1) +
  scale_color_manual(values = c(
    "TPL/DRB.control" = "red",
    "TPL/DRB.treatment" = "blue",
    "U1 AMO.control" = "red",
    "U1 AMO.treatment" = "blue"
  )) +
  ylim(c(10, 110)) +
  theme_minimal(base_size = 20) +
  labs(x = NULL) +
  ylab(bquote('Volume '(µm^3))) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 2.3, lineend = "round"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line = element_line(colour = "black", linewidth = 2.3, lineend = "round"),
    axis.text = element_text(size = 40, colour = "black"),
    axis.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 30)
  ) +
  theme(legend.text = element_text(size = 25),
        legend.title = element_text(size = 30)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  
  theme(legend.position = "none") 
dev.off()


# Filter for control samples only
control_data <- w %>%
  filter(sample == "control")

# Get the two treatment groups
treat1_data <- control_data$Freq[control_data$treatment == "TPL/DRB"]
treat2_data <- control_data$Freq[control_data$treatment == "U1 AMO"]  

# Perform KS test
ks_result <- ks.test(treat1_data, treat2_data)

# Print results
cat("KS test between control samples:\n")
cat("TPL/DRB vs", unique(control_data$treat_id)[2], "\n")
cat("D statistic:", ks_result$statistic, "\n")
cat("p-value:", ks_result$p.value, "\n")

# Perform KS tests for each treatment
ks_results <- w %>%
  group_by(treatment) %>%
  summarise(
    ks_statistic = ks.test(Freq[sample == "control"], 
                           Freq[sample == "treatment"])$statistic,
    ks_pvalue = ks.test(Freq[sample == "control"], 
                        Freq[sample == "treatment"])$p.value,
    .groups = "drop"
  )

print(ks_results)



###################   Plotting Split Violins  ###############


##########  Load data

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


rm(DGE_HS42_V1_intronTCReads,DGE_intronTCReads,DGE_HS42_V1_exonTCReads, DGE_exonTCReads,DGE_HS42_V1_allTCReads,DGE_allTCReads)






#### Fig. 2d - Allreads

###### Define tables to plot
treatments_keep <- c("Transcription_block_allReads",
                     "U1AMO_4h_30uM_allReads")

combined_sub <- DGE_allTCReads_combined %>% filter(treatment %in% treatments_keep)
sig_only_sub <- combined_sub %>% filter(Significance == "Significant")

treatment_order <- c("Transcription_block_allReads",
                     "U1AMO_4h_30uM_allReads")  

combined_sub$treatment <- factor(combined_sub$treatment, levels = treatment_order)
sig_only_sub$treatment <- factor(sig_only_sub$treatment, levels = treatment_order)

treatment_levels <- levels(combined_sub$treatment)

treatment_levels <- levels(factor(combined_sub$treatment))

## Helper functions
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

compute_boxplot_stats <- function(df, side, group_col = "treatment") {
  treatments <- unique(df[[group_col]])
  result <- lapply(treatments, function(trt) {
    vals <- df$Log2FC[df[[group_col]] == trt]
    vals <- vals[!is.na(vals)]
    x_center <- as.numeric(factor(trt, levels = treatment_levels))
    x_offset <- if (side == "left") -0.1 else 0.1
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
      y     = 4,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  })
  bind_rows(result)
}

### Compute plot data
scale_factor <- 0.4

left_dens  <- compute_half_violin(combined_sub, "left")
right_dens <- compute_half_violin(sig_only_sub, "right")

left_dens$x_pos  <- as.numeric(factor(left_dens$treatment,  levels = treatment_levels)) + left_dens$density  * scale_factor
right_dens$x_pos <- as.numeric(factor(right_dens$treatment, levels = treatment_levels)) + right_dens$density * scale_factor

box_left  <- compute_boxplot_stats(combined_sub, "left")
box_right <- compute_boxplot_stats(sig_only_sub, "right")
box_width <- 0.04

annot_df <- compute_annotations(combined_sub)

### Plot
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
pdf("TxBlock_vs_U1AMO_SplitViolin_SLAMseq_allreads_v10.pdf", height = 6, width = 6.7)
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
  scale_x_continuous(breaks = seq_along(treatment_levels), labels = treatment_levels) +
  scale_fill_manual(values = c("All" = "cornflowerblue", "Significant" = "coral1")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.ticks.x  = element_line(linewidth = 2.5, lineend = "round"), 
    axis.ticks.y  = element_line(linewidth = 2.5, lineend = "round"), 
    axis.ticks.length = unit(0.3, "cm"),
    axis.line     = element_line(linewidth = 2.5, lineend = "round"),
    axis.text.y = element_text(size = 46, margin = margin(r = 10)),
    axis.title.y  = element_blank(),   
    legend.position = "none"
  ) +
  labs(y = "Log2 Fold Change")
dev.off()








#####
#### Fig. 2e - introns
###### Define tables to plot

treatments_keep <- c("Transcription_block_intronsOnly",
                     "U1AMO_4h_30uM_intronsOnly")

combined_sub <- DGE_intronTCReads_combined %>% filter(treatment %in% treatments_keep)
sig_only_sub <- combined_sub %>% filter(Significance == "Significant")
#
treatment_order <- c("Transcription_block_intronsOnly",
                     "U1AMO_4h_30uM_intronsOnly")  
combined_sub$treatment <- factor(combined_sub$treatment, levels = treatment_order)
sig_only_sub$treatment <- factor(sig_only_sub$treatment, levels = treatment_order)

treatment_levels <- levels(combined_sub$treatment)

treatment_levels <- levels(factor(combined_sub$treatment))


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
      y     = 7,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  })
  bind_rows(result)
}

### Compute plot data
scale_factor <- 0.4

left_dens  <- compute_half_violin(combined_sub, "left")
right_dens <- compute_half_violin(sig_only_sub, "right")

left_dens$x_pos  <- as.numeric(factor(left_dens$treatment,  levels = treatment_levels)) + left_dens$density  * scale_factor
right_dens$x_pos <- as.numeric(factor(right_dens$treatment, levels = treatment_levels)) + right_dens$density * scale_factor

box_left  <- compute_boxplot_stats(combined_sub, "left")
box_right <- compute_boxplot_stats(sig_only_sub, "right")
box_width <- 0.04

annot_df <- compute_annotations(combined_sub)

### Plot
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
pdf("TxBlock_vs_U1AMO_SplitViolin_SLAMseq_introns.pdf", height = 6, width = 6.7)
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
  scale_x_continuous(breaks = seq_along(treatment_levels), labels = treatment_levels) +
  scale_fill_manual(values = c("All" = "cornflowerblue", "Significant" = "coral1")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme_classic(base_size = 14) + 
  theme(
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.ticks.x  = element_line(linewidth = 2.5, lineend = "round"),
    axis.ticks.y  = element_line(linewidth = 2.5, lineend = "round"), 
    axis.ticks.length = unit(0.3, "cm"),
    axis.line     = element_line(linewidth = 2.5, lineend = "round"),
    axis.text.y = element_text(size = 46, margin = margin(r = 10)),
    axis.title.y  = element_blank(), 
    legend.position = "none"
  ) +
  labs(y = "Log2 Fold Change")
dev.off()




#### ED Fig. 2b - exons

treatments_keep <- c("Transcription_block_exonsOnly",
                     "U1AMO_4h_30uM_exonsOnly")

combined_sub <- DGE_exonTCReads_combined %>% filter(treatment %in% treatments_keep)
sig_only_sub <- combined_sub %>% filter(Significance == "Significant")

treatment_order <- c("Transcription_block_exonsOnly",
                     "U1AMO_4h_30uM_exonsOnly")
combined_sub$treatment <- factor(combined_sub$treatment, levels = treatment_order)
sig_only_sub$treatment <- factor(sig_only_sub$treatment, levels = treatment_order)


treatment_levels <- levels(combined_sub$treatment)

treatment_levels <- levels(factor(combined_sub$treatment))


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
      y     = 7,
      label = paste0("total=", n_total, "\nup=", n_up, "\ndown=", n_down)
    )
  })
  bind_rows(result)
}

### Compute plot data
scale_factor <- 0.4

left_dens  <- compute_half_violin(combined_sub, "left")
right_dens <- compute_half_violin(sig_only_sub, "right")

left_dens$x_pos  <- as.numeric(factor(left_dens$treatment,  levels = treatment_levels)) + left_dens$density  * scale_factor
right_dens$x_pos <- as.numeric(factor(right_dens$treatment, levels = treatment_levels)) + right_dens$density * scale_factor

box_left  <- compute_boxplot_stats(combined_sub, "left")
box_right <- compute_boxplot_stats(sig_only_sub, "right")
box_width <- 0.04

annot_df <- compute_annotations(combined_sub)

### Plot
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
pdf("TxBlock_vs_U1AMO_SplitViolin_SLAMseq_exonReads_v10.pdf", height = 6, width = 6.7)
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
  scale_x_continuous(breaks = seq_along(treatment_levels), labels = treatment_levels) +
  scale_fill_manual(values = c("All" = "cornflowerblue", "Significant" = "coral1")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  coord_cartesian(ylim = c(-7, 7)) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x   = element_blank(),
    axis.title.x  = element_blank(),
    axis.ticks.x  = element_line(linewidth = 2.5, lineend = "round"), 
    axis.ticks.y  = element_line(linewidth = 2.5, lineend = "round"),  
    axis.ticks.length = unit(0.3, "cm"),
    axis.line     = element_line(linewidth = 2.5, lineend = "round"),  
    axis.text.y = element_text(size = 46, margin = margin(r = 10)),
    axis.title.y  = element_blank(),  
    legend.position = "none"
  ) +
  labs(y = "Log2 Fold Change")
dev.off()


####### Figure 2f-h 
######### Crane maps were generated by "make_cranemaps_Fig2.ipynb"

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
seg_size <- 30
u <- 2.4
lw <- 4.2

## Fig.2d - PC1-colored track (chr11)
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



## Fig.2f - LOS residuals for NT1 chr 11
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")

x <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/150kb/LOS_residuals/LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103.bedGraph",
                sep="\t", header=TRUE)
x <-x[x$chrom == "chr11", ]
pdf("LOSres_NT1_track_ch11_xaxis_150kb.pdf", height = 2.4, width = 9)
par(mfrow=c(1,1), mar=c(2, 4, 2, 2) + 0.1)

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=FALSE) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()

#####

## Fig.2g - LOS residuals for TD1 chr 11
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
y <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/150kb/LOS_residuals/LOS_residuals_range2Mb_K562_PIIBTD240mPD_DpnII_R1_20200103.bedGraph",
                sep="\t", header=TRUE)
y <-y[y$chrom == "chr11", ]
pdf("LOSres_TD1_track_ch11_xaxis_150kb.pdf", height = 2.4, width = 9)
par(mfrow=c(1,1), mar=c(2, 4, 2, 2) + 0.1)

plot(y$start/1000000, y$LOS_residuals_range2Mb_K562_PIIBTD240mPD_DpnII_R1_20200103, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=FALSE) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()

#####



## Fig.2h - Difference in LOS residuals betwee NT1 and TD1 chr 11
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
z <- merge(x,y, by = c("chrom","start","end"),all.x= TRUE)
z$dif <- z$LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103 - z$LOS_residuals_range2Mb_K562_PIIBTD240mPD_DpnII_R1_20200103
z$chrom <- factor(z$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
z <- z[order(z$chrom, z$start,z$end),]
rownames(z) <- NULL

pdf("LOSres_NTTD1dif_track_ch11_xaxis_150kb.pdf", height = 2.5, width = 9)
par(mfrow=c(1,1), mar=c(2.5, 4, 2, 2) + 0.1)

plot(z$start/1000000, z$dif, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=TRUE,mgp=c(3, 1.5, 0)) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )


dev.off()


#####





## Fig.2i - zoom chr3 for TPL/DRB and U1 treatments

seg_size <- 25

begin <- 33000000
term <- 60000000

t <- 1.3
lw  <- 2.7

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")
x <- df_residuals[df_residuals$chrom == "chr3",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_TREATdif_track_ch3_zoom_TxU14hRNase_wideA2_alt.pdf", height = 2.7, width = 4.7)
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
lw  <- 2.6
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/")

x <- df_residuals[df_residuals$chrom == "chr3",]
x <- x[x$start>begin & x$end<term,]
pdf("PC1_SPINsegs_TREATdif_track_ch3_zoom_TxU14h_wideB2_v10.pdf", height = 2.85, width = 4.7)
par(mfrow=c(2,1), mar=c(0.75, 4, 1.8, 2) + 0.1,oma=c(1, 0, 0, 0), mgp=c(3, 0.6, 0))
plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif, type="l", col="black",
     ylim=c(-0.22, 0.15), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=t, tck=-0.08,labels=FALSE)
axis(2, lwd=lw , cex.axis=t, tck=-0.08, mgp=c(3, 0.6, 0))
box(bty="l", lwd=lw )

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif, type="l", col="black",
     ylim=c(-0.22, 0.15), xlab = "", ylab= "", axes=FALSE, lwd=1)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=t, tck=-0.08,labels=TRUE)
axis(2, lwd=lw , cex.axis=t, tck=-0.08, mgp=c(3, 0.5, 0))
box(bty="l", lwd=lw )

dev.off()




## Fig.2i - Global violins 

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/violin")


NT <- df_residuals[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif','LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif','LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628_dif')]


NTTD <- NT[,c(4, 6)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "TPL/DRB"
NTTD <- NTTD[NTTD$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN2_global_NTTD1_violin_alt2.pdf", width=5, height=2.1)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#ABDDA4",  "#5E4FA2")) +
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


aggregate(LOS_residuals ~ SPIN, data = na.omit(NTTD), FUN = mean)


U14h <- NT[,c(4, 7)]
colnames(U14h)[2] <- "LOS_residuals"
U14h$treat <- "U1 AMO"
U14h <- U14h[U14h$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNAtreat/violin")
pdf("LOSr_selectSPIN_global_U14h_violin_alt2.pdf", width=5, height=2.1)
ggplot(na.omit(U14h), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142", "#F46D43","#FDAE61", "#ABDDA4",  "#5E4FA2")) +
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

aggregate(LOS_residuals ~ SPIN, data = na.omit(U14h), FUN = mean)


