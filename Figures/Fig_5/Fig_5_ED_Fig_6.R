library(tidyverse)
library(rtracklayer)
library(data.table)

options(scipen = 999)


df_residuals <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/50kb/dataframes/df_residuals_20260812_50kb.bed",
                           sep="\t", header=TRUE)
df_residuals$seg_1 <- -1
df_residuals <- cbind(df_residuals, replicate(8, df_residuals$seg_1))
n <- ncol(df_residuals)
names(df_residuals)[(n-7):n] <- c("Speckle_seg", "Interior_Act1_seg", "Interior_Act2_seg", "Interior_Repr1_seg","Interior_Repr2_seg","Near_Lm1_seg","Near_Lm2_seg","Lamina_seg")

df_residuals <- df_residuals %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),
                  ordered = TRUE))


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


df_residuals_150kb_gc <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/150kb/dataframes/df_residuals_20260326_150kb_GC_alt.bed",
                                    sep="\t", header=TRUE)
df_residuals_150kb_gc$chrom <- factor(df_residuals_150kb_gc$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals_150kb_gc <- df_residuals_150kb_gc[order(df_residuals_150kb_gc$chrom, df_residuals_150kb_gc$start,df_residuals_150kb_gc$end),]

df_residuals_150kb_gc <- df_residuals_150kb_gc %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),
                  ordered = TRUE))

df_residuals_150kb_gc["LOS_residuals_range2Mb_K562_RI2NT1_DpnII_R1_20200103_dif"] <- df_residuals_150kb_gc["LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103"] - df_residuals_150kb_gc["LOS_residuals_range2Mb_K562_RNaseIN4hPD_DpnII_R2_20230628"]
df_residuals_150kb_gc["LOS_residuals_range2Mb_K562_NT1RA2_DpnII_R1_20200103_dif"] <- df_residuals_150kb_gc["LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103"] - df_residuals_150kb_gc["LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628"]



######### df_residuals at 1Mb
df_residuals_1000kb_gc <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/1000kb/dataframes/df_residuals_20260326_1000kb_GC_alt.bed",
                                     sep="\t", header=TRUE)



df_residuals_1000kb_gc$chrom <- factor(df_residuals_1000kb_gc$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
df_residuals_1000kb_gc <- df_residuals_1000kb_gc[order(df_residuals_1000kb_gc$chrom, df_residuals_1000kb_gc$start,df_residuals_1000kb_gc$end),]

df_residuals_1000kb_gc <- df_residuals_1000kb_gc %>%
  mutate(SPIN = recode(SPIN, "Interior_Act3" = "Interior_Act2") %>%
           factor(levels = c("Speckle", "Interior_Act1", "Interior_Act2",
                             "Interior_Repr1", "Interior_Repr2",
                             "Near_Lm1", "Near_Lm2", "Lamina"),
                  ordered = TRUE))




SONtsaseq <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/DL/tsaseq_1mb.bedgraph",
                        sep="\t", header=FALSE)
colnames(SONtsaseq) <- c("chrom","start","end","SON_TSAseq_signal")

df_residuals_1000kb_gc  <- merge(df_residuals_1000kb_gc, SONtsaseq, by = c("chrom","start","end"), all.x = TRUE)

df_residuals_1000kb_gc["LOS_residuals_range2Mb_K562_NT1RI2_DpnII_R1_20200103_dif"] <- df_residuals_1000kb_gc["LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103"] - df_residuals_1000kb_gc["LOS_residuals_range2Mb_K562_RNaseIN4hPD_DpnII_R2_20230628"]
df_residuals_1000kb_gc["LOS_residuals_range2Mb_K562_NT1RA2_DpnII_R1_20200103_dif"] <- df_residuals_1000kb_gc["LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103"] - df_residuals_1000kb_gc["LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628"]


##########  Load DGE data

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



###### Fig. 5a-c

######################  Correlate LOS destabilization with differential expression   ####################

###################### Calculate number of genes up/down-regulated after each treatment (all genes)  ############

gtf <- import("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/ENSEMBL_genes/Homo_sapiens.GRCh38.113.gtf.gz")
genes_gtf <- as.data.frame(gtf[gtf$type == "gene"])


##### Prepare TSS from your GTF table 
tss <- genes_gtf %>%
  filter(type == "gene") %>%
  mutate(
    chr    = paste0("chr", seqnames),   
    TSS    = ifelse(strand == "+", start, end)
  ) %>%
  select(gene_id, gene_name, chr, TSS, strand)

#### Load your SPIN bed file
SPIN <- read.table("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Sequencing/data/K562_data/SPIN/hg38/SPIN_50kb_rebinned.bed", sep="\t",
                   header=TRUE)
SPIN <- SPIN[, !names(SPIN) %in% c("n_bins", "n_unique")]
SPIN$SPIN <- factor(SPIN$SPIN, order = TRUE, levels =c("Speckle", "Interior_Act1", "Interior_Act2", "Interior_Act3", "Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina"))
#### Annotate each TSS with SPIN state
# data.table non-equi join: find the bin where TSS falls within [start, end]
setDT(tss)
setDT(SPIN)

tss_spin <- SPIN[tss,
                 on = .(chrom = chr, start <= TSS, end >= TSS),
                 .(gene_id, SPIN),
                 nomatch = NA
]

### Make lookup table and annotate DGE tables
spin_lookup <- tss_spin %>%
  select(gene_id, SPIN) %>%
  filter(!is.na(SPIN)) 




DGE_exonTCReads_combined_all <- DGE_exonTCReads_combined %>%
  left_join(spin_lookup, by = c("Gene" = "gene_id"))

DGE_exonTCReads_combined_all <- DGE_exonTCReads_combined_all %>%
  left_join(tss, by = c("Gene" = "gene_id"))


DGE_allTCReads_combined_all <- DGE_allTCReads_combined %>%
  left_join(spin_lookup, by = c("Gene" = "gene_id"))

DGE_allTCReads_combined_all <- DGE_allTCReads_combined_all %>%
  left_join(tss, by = c("Gene" = "gene_id"))

##### Speckles only

DGE_exonTCReads_combined_sp <- DGE_exonTCReads_combined_all[DGE_exonTCReads_combined_all$SPIN == "Speckle",]

DGE_exonTCReads_combined_NOsp <- DGE_exonTCReads_combined_all[DGE_exonTCReads_combined_all$SPIN %in% c("Interior_Act1", "Interior_Act2", "Interior_Act3"),]

DGE_exonTCReads_combined_Repr <- DGE_exonTCReads_combined_all[DGE_exonTCReads_combined_all$SPIN %in% c("Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina"), ]


treatments_keep <- c(
  "Transcription_block_exonsOnly",
  "U1AMO_4h_30uM_exonsOnly",
  "Heat_shock_80min_exonsOnly",
  "Speckle-dTAGvsNEG_6h_exonsOnly"
)

name_map <- c(
  "Transcription_block_exonsOnly"   = "tx_block",
  "U1AMO_4h_30uM_exonsOnly"         = "U1_AMO",
  "Heat_shock_80min_exonsOnly"      = "heat_shock",
  "Speckle-dTAGvsNEG_6h_exonsOnly"  = "dTAGvsNeg"
)


count_down_by_treatment <- function(df, treatments_keep, name_map) {
  df %>%
    filter(treatment %in% names(name_map)) %>%
    mutate(treatment = factor(treatment, levels = names(name_map))) %>%
    group_by(treatment, .drop = FALSE) %>%
    summarise(n_down = sum(Color_factor == "Significantly_Down", na.rm = TRUE),
              .groups = "drop") %>%
    mutate(treatment = factor(name_map[as.character(treatment)], levels = name_map))
}

DGE_count_exons <- count_down_by_treatment(
  DGE_exonTCReads_combined, treatments_keep, name_map
)

DGE_count_exons_NOsp <- count_down_by_treatment(
  DGE_exonTCReads_combined_NOsp, treatments_keep, name_map
)

DGE_count_exons_sp <- count_down_by_treatment(
  DGE_exonTCReads_combined_sp, treatments_keep, name_map
)

DGE_count_exons_Repr  <- count_down_by_treatment(
  DGE_exonTCReads_combined_Repr, treatments_keep, name_map
)



#####

############## Calculating the avg diff in LOS residuals at speckles for each  treatment ####

x <- df_residuals %>%
  filter(SPIN == "Speckle") %>%
  select(
    LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif,
    LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif,
    LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif,
    LOS_residuals_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif
  )

LOSres_avg_SpSPIN <- x %>%
  summarise(across(where(is.numeric), \(v) mean(v, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "mean_LOS_residual")

############## Calculating the avg diff in LOS residuals at non-speckle active regions for each  treatment ####

x <- df_residuals %>%
  filter(SPIN %in% c("Interior_Act1", "Interior_Act2", "Interior_Act3")) %>%
  select(
    LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif,
    LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif,
    LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif,
    LOS_residuals_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif
  )

LOSres_avg_ActSPIN <- x %>%
  summarise(across(where(is.numeric), \(v) mean(v, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "mean_LOS_residual")


############## Calculating the avg diff in LOS residuals at inactive regions for each  treatment ####

x <- df_residuals %>%
  filter(SPIN %in% c("Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina")) %>%
  select(
    LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif,
    LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif,
    LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif,
    LOS_residuals_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif
  )

LOSres_avg_ReprSPIN <- x %>%
  summarise(across(where(is.numeric), \(v) mean(v, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "sample", values_to = "mean_LOS_residual")



###############################  PLotting LOS_diff vs DGE scatters #########################

### Fig.5a

y <- cbind(DGE_count_exons,LOSres_avg_SpSPIN)
rownames(y) <- NULL
y$treat_name <- c("Tx block","U1 AMO", "Heat shock","Sp-dTAG")

#### PLot scatter for all DE genes
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/scatter")
pdf("LOS_SLAM_correlation_scatter_exonReads.pdf", height = 3, width = 3.5)
ggplot(y, aes(x = n_down, y = mean_LOS_residual, label = treat_name)) +
  geom_point(size=3) +
  geom_text(vjust = -1, hjust = 0.35, color="black", size=5) + 
  labs(x = "Down-regulated genes",
       y = "Mean difference in LOS residuals at Speckles") +
  xlim(-1500, 11000) +
  ylim(-0.046, 0.01) +
  theme_classic() +
  theme(axis.ticks = element_line(colour = "black", linewidth = 1.3, lineend = "round")) +
  theme(axis.ticks.length = unit(0.25, "cm")) +
  theme(axis.line = element_line(colour = "black", linewidth = 1.3, lineend = "round")) +
  theme(axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10)), 
        axis.text.y = element_text(size = 22, colour = "black"),
        axis.title = element_text(size = 30, face = "bold")) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) 


dev.off()




### Fig.5b

y <- cbind(DGE_count_exons_sp,LOSres_avg_SpSPIN)
rownames(y) <- NULL
y$treat_name <- c("Tx block","U1 AMO", "Heat shock","Sp-dTAG")

#### PLot scatter for only DE genes within speckle-associated regions
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/scatter")
pdf("LOS_SLAM_correlation_scatter_exonReads_sp.pdf", height = 3, width = 2.7)
ggplot(y, aes(x = n_down, y = mean_LOS_residual, label = treat_name)) +
  geom_point(size=3) +
  geom_text(vjust = -1, hjust = 0.35, color="black", size=5) + 
  labs(x = "Down-regulated genes",
       y = "Mean difference in LOS residuals at Speckles") +
  scale_x_continuous(
    breaks = seq(0, 2000, by = 1000),
    limits = c(-350, 2550)
  ) +
  ylim(-0.046, 0.01) +
  theme_classic() +
  theme(axis.ticks = element_line(colour = "black", linewidth = 1.3, lineend = "round"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.line = element_line(colour = "black", linewidth = 1.3, lineend = "round"),
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10)),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


dev.off()



### Fig.5c

y <- cbind(DGE_count_exons_NOsp,LOSres_avg_ActSPIN)
rownames(y) <- NULL
y$treat_name <- c("Tx block","U1 AMO", "Heat shock","Sp-dTAG")

#### PLot scatter for DE genes within Int Active regions
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/scatter")
pdf("LOS_SLAM_correlation_scatter_exonReads_NOsp.pdf", height = 3, width = 2.7)
ggplot(y, aes(x = n_down, y = mean_LOS_residual)) +
  geom_point(size=3) +

  labs(x = "Down-regulated genes",
       y = "Mean difference in LOS residuals at Speckles") +
  scale_x_continuous(
    breaks = seq(0, 6000, by = 3000),
    limits = c(-800, 6000)
  ) +
  ylim(-0.046, 0.01) +
  theme_classic() +
  theme(axis.ticks = element_line(colour = "black", linewidth = 1.3, lineend = "round"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.line = element_line(colour = "black", linewidth = 1.3, lineend = "round"),
        axis.text.x = element_text(size = 22, colour = "black", margin = margin(t = 10)),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 30, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


dev.off()


####### Fig. 5e 
######### Crane maps in Fig. 5e and Fig. 5g were generated using 'make_cranemaps_Fig5.ipynb'

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")
seg_size <- 30
u <- 2.4
lw <- 4.2


## PC1-colored track (chr3)
x <-df_residuals[df_residuals$chrom == "chr3", ]
pdf("PC1_SPIN_segs_track_ch3_xaxis.pdf", height = 4.2, width = 9)
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

## Fig.5e- LOS residuals for NT1 chr3
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

x <- df_residuals_150kb_gc
x <-x[x$chrom == "chr3", ]
pdf("LOSres_NT1_track_ch3_xaxis_150kb.pdf", height = 2.5, width = 9)
par(mfrow=c(1,1), mar=c(2.5, 4, 2, 2) + 0.1,mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=FALSE,mgp=c(3, 1.5, 0)) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()

#####

## Fig.5e- LOS residuals for RA2 chr 11
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

x <- df_residuals_150kb_gc
x <-x[x$chrom == "chr3", ]
pdf("LOSres_RA2_track_ch11_xaxis_150kb.pdf", height = 2.5, width = 9)
par(mfrow=c(1,1), mar=c(2.5, 4, 2, 2) + 0.1,mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628, type="l", col="black",
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=FALSE,mgp=c(3, 1.5, 0))
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()




## Fig.5e- - Difference in LOS residuals (NT1RA2) for chr 3
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

x <- df_residuals_150kb_gc
x <-x[x$chrom == "chr3", ]
pdf("LOSres_NT1RA2diff_track_ch3_xaxis_150kbx.pdf", height = 2.5, width = 9)
par(mfrow=c(1,1), mar=c(2.5, 4, 2, 2) + 0.1,mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_NT1RA2_DpnII_R1_20200103_dif, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=FALSE,mgp=c(3, 1.5, 0)) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()



################## Violins for specific SPIN states (Speckle, Interior_Act1, Interio_repr2, Lamina) ###############################################

NT <- df_residuals_150kb_gc[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif',
                               'LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif','LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628_dif',
                               'LOS_residuals_range2Mb_K562_RI2NT1_DpnII_R1_20200103_dif','LOS_residuals_range2Mb_K562_NT1RA2_DpnII_R1_20200103_dif')]


### Fig.5f

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

####### NT1 - RA2 
NTTD <- NT[,c(4, 10)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "NTvRI"
NTTD <- NTTD[NTTD$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN_global_NT1RA2_violin_alt2.pdf", width=6, height=3.05)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61","#ABDDA4",  "#5E4FA2")) +
  ylim(-0.12, 0.1) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 1.5, lineend = "round"),
        axis.ticks = element_line(linewidth = 1.5, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"),  
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.y = element_text(size = 28, colour = "black", margin = margin(r = 5, unit = "pt")))

dev.off()


options(scipen = 999)
aggregate(LOS_residuals ~ SPIN, data = na.omit(NTTD), FUN = mean)



## Fig.5g - SPIN seg

## (chr3)
x <-df_residuals[df_residuals$chrom == "chr3", ]
pdf("SPIN_segs_track_ch3_xaxis.pdf", height = 4.2, width = 9)
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
axis(1, lwd=lw, cex.axis=2.2, labels=FALSE) 


dev.off()



## Fig.5g - LOS residuals for RI2 chr3
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

x <- df_residuals_150kb_gc
x <-x[x$chrom == "chr3", ]
pdf("LOSres_RI2_track_ch3_xaxis_150kb.pdf", height = 2.4, width = 9)
par(mfrow=c(1,1), mar=c(2, 4, 2, 2) + 0.1,mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_RNaseIN4hPD_DpnII_R2_20230628, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=FALSE) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()



## Fig.5g - Difference in LOS residuals for RI2NT1 chr 3
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

x <- df_residuals_150kb_gc
x <-x[x$chrom == "chr3", ]
pdf("LOSres_RI2NT1diff_track_ch3_xaxis_150kbx.pdf", height = 2.55, width = 9)
par(mfrow=c(1,1), mar=c(2.5, 4, 2, 2) + 0.1,mgp=c(3, 1.5, 0))

plot(x$start/1000000, x$LOS_residuals_range2Mb_K562_RI2NT1_DpnII_R1_20200103_dif, type="l", col="black", 
     ylim=c(-0.18, 0.1), xlab = "", ylab= "", axes=FALSE, lwd=2)
abline(h = 0,col = "black", lty=2,lwd=2)
axis(1, lwd=lw , cex.axis=2.2, labels=TRUE,mgp=c(3, 1.5, 0)) 
axis(2, lwd=lw , cex.axis=u)
box(bty="l", lwd=lw )

dev.off()



### Fig.5h

setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

####### NT1 - RI2 
NTTD <- NT[,c(4, 9)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "NT1vsRA2"
NTTD <- NTTD[NTTD$SPIN == c("Speckle", "Interior_Act1","Interior_Act2","Interior_Repr2","Lamina"),]

pdf("LOSr_selectSPIN_global_RI2NT1_violin_alt2.pdf", width=4.9, height=3.05)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142",  "#F46D43","#FDAE61", "#ABDDA4",  "#5E4FA2")) +
  ylim(-0.12, 0.1) +
  geom_hline(yintercept = 0,color="black",linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 1.5, lineend = "round"),
        axis.ticks = element_line(linewidth = 1.5, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"),  
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

dev.off()

options(scipen = 999)
aggregate(LOS_residuals ~ SPIN, data = na.omit(NTTD), FUN = mean)



##### Extended Data FIg. 6a-d

################################################################################################################
###################### Colored LOS vs SLAM CPM scatters ########################################################


NT <- df_residuals_1000kb_gc
NT_spin <- NT[NT$SPIN %in% c("Speckle", "Interior_Act1"), ]

#### Extended Data Fig. 6a

#######################  Difference in LOS residuals on the y axis (after tx block) ############

x_var <- colnames(NT)[130]
y_var <- colnames(NT)[105]
c_var <- colnames(NT)[5]

NT <- NT[!is.na(NT[[x_var]]) & !is.na(NT[[y_var]]), ]



setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

pdf("LOSdif_vs_SLAMCPSp_SPINColor_Speckle_NTTD1.pdf", height = 4, width = 4)
ggplot() +
  geom_point(data = NT,
             aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[c_var]]),
             size = 2) +
  scale_color_manual(values = c("Speckle" = "#9E0142", "Interior_Act1" ="#F46D43", "Interior_Act2"="#FDAE61", "Interior_Repr1"="#FFFFBF", "Interior_Repr2"="#ABDDA4", 
                                "Near_Lm1"="#66C2A5", "Near_Lm2"="#3288BD", "Lamina"="#5E4FA2")) +
  scale_x_log10(limits = c(100, 1e6)) +
  geom_vline(xintercept = (quantile(NT$HS37_20241101_1000000_SpNorm, 0.85, na.rm = TRUE)), color = "black", linetype = "dashed", linewidth = 0.8) +
  labs(x = x_var, y = y_var, title = NULL) +
  theme_classic() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_text(size = 13),
    axis.line       = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks      = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.margin     = margin(t = 5.5, r = 20, b = 5.5, l = 5.5)
  )
dev.off()



###### Extended Data Fig. 6c


###########  Plot LOSdif NTTD1 vs TSASON  - color by SPIN (all SPIN)

####### Difference in LOS residuals (NTTD1) ########
x_var <- colnames(NT)[154]
y_var <- colnames(NT)[105]
c_var <- colnames(NT)[5]

NT <- NT[!is.na(NT[[x_var]]) & !is.na(NT[[y_var]]), ]


# Global GC range across all SPIN states
gc_min <- min(NT[[c_var]], na.rm = TRUE)
gc_max <- max(NT[[c_var]], na.rm = TRUE)
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

pdf("LOSdiff_vs_TSASON_SPINColor_allSPIN_NTTD1.pdf", height = 4, width = 4)
ggplot() +
  geom_point(data = NT,
             aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[c_var]]),
             size = 2) +
  scale_color_manual(values = c("Speckle" = "#9E0142", "Interior_Act1" ="#F46D43", "Interior_Act2"="#FDAE61", "Interior_Repr1"="#FFFFBF", "Interior_Repr2"="#ABDDA4", 
                                "Near_Lm1"="#66C2A5", "Near_Lm2"="#3288BD", "Lamina"="#5E4FA2")) +
  labs(x = x_var, y = y_var, title = NULL) +
  scale_x_continuous(breaks = seq(-1, 3, by = 1)) +
  scale_y_continuous(breaks = seq(-0.2, 0.1, by = 0.1)) +
  xlim(-1.5,3) +
  theme_classic() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_text(size = 13),
    axis.line       = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks      = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.margin     = margin(t = 5.5, r = 20, b = 5.5, l = 5.5)
  )
dev.off()



######### Extended Data Fig. 6d

###########  Plot GCcontent vs TSASON  - color by SPIN (all SPIN)

####### Difference in LOS residuals (NT1RA2) ########
x_var <- colnames(NT)[154]
y_var <- colnames(NT)[153]
c_var <- colnames(NT)[5]

NT <- NT[!is.na(NT[[x_var]]) & !is.na(NT[[y_var]]), ]

gc_min <- min(NT[[c_var]], na.rm = TRUE)
gc_max <- max(NT[[c_var]], na.rm = TRUE)
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

pdf("GCcontent_vs_TSASON_SPINColor_allSPIN.pdf", height = 4, width = 4)
ggplot() +
  geom_point(data = NT,
             aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[c_var]]),
             size = 2) +
  scale_color_manual(values = c("Speckle" = "#9E0142", "Interior_Act1" ="#F46D43", "Interior_Act2"="#FDAE61", "Interior_Repr1"="#FFFFBF", "Interior_Repr2"="#ABDDA4", 
                                "Near_Lm1"="#66C2A5", "Near_Lm2"="#3288BD", "Lamina"="#5E4FA2")) +
  labs(x = x_var, y = y_var, title = NULL) +
  scale_x_continuous(breaks = seq(-1, 3, by = 1)) +
  xlim(-1.5,3) +
  theme_classic() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_text(size = 13),
    axis.line       = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks      = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.margin     = margin(t = 5.5, r = 20, b = 5.5, l = 5.5)
  )
dev.off()


##### Extended Data Fig. 6b

################################# Subset quartile of highest expression in NT table and plot  #######################  


NT <- df_residuals_1000kb_gc
NT$quantile_CPM <- ntile(NT$HS37_20241101_1000000_SpNorm, 100)
NT_highCPM <- NT[NT$quantile_CPM >= 85,]

####### Difference in LOS residuals (NT1RA2) ########
x_var <- colnames(NT)[154]
y_var <- colnames(NT)[105]
c_var <- colnames(NT)[5]

NT <- NT[!is.na(NT[[x_var]]) & !is.na(NT[[y_var]]), ]

gc_min <- min(NT[[c_var]], na.rm = TRUE)
gc_max <- max(NT[[c_var]], na.rm = TRUE)
setwd("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/Figures/Figure_RNA/")

pdf("LOSdiff_vs_TSASON_SPINColor_highCPM_NT1RI2.pdf", height = 4, width = 4)
ggplot() +
  geom_point(data = NT_highCPM,
             aes(x = .data[[x_var]], y = .data[[y_var]], color = .data[[c_var]]),
             size = 2) +
  scale_color_manual(values = c("Speckle" = "#9E0142", "Interior_Act1" ="#F46D43", "Interior_Act2"="#FDAE61", "Interior_Repr1"="#FFFFBF", "Interior_Repr2"="#ABDDA4", 
                                "Near_Lm1"="#66C2A5", "Near_Lm2"="#3288BD", "Lamina"="#5E4FA2")) +
  labs(x = x_var, y = y_var, title = NULL) +
  scale_x_continuous(breaks = seq(-1, 3, by = 1)) +
    xlim(-1.5,3) +
  theme_classic() +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_text(size = 13),
    axis.line       = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks      = element_line(linewidth = 0.8, lineend = "round"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    plot.title      = element_text(size = 12, face = "bold"),
    plot.margin     = margin(t = 5.5, r = 20, b = 5.5, l = 5.5)
  )
dev.off()




################ Glocal Violins ##########################


NT <- df_residuals_150kb_gc[,c('chrom' , 'start', 'end', 'SPIN','PC1','LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif',
                               'LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif','LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628_dif',
                               'LOS_residuals_range2Mb_K562_RI2NT1_DpnII_R1_20200103_dif','LOS_residuals_range2Mb_K562_NT1RA2_DpnII_R1_20200103_dif')]

### Extended Data Fig. 6g

####### NT1 - RA2 
NTTD <- NT[,c(4, 10)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "NT1vsRA2"

#Violin plotting the LOS residuals for dTAG 6h
pdf("LOSr_allSPIN_global_NT1RA2_violin_alt2.pdf", width=10, height=3)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.215, 0.145) +
  geom_hline(yintercept = 0, color="black", linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2, lineend = "round"),
        axis.ticks = element_line(linewidth = 2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm")) +
  theme(axis.title.y=element_blank())

dev.off()


### Extended Data Fig. 6h

####### RI2 - NT1 
NTTD <- NT[,c(4, 9)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "RIvsNT"

#Violin plotting the LOS residuals for dTAG 6h
pdf("LOSr_allSPIN_global_RI2NT1_violin_alt2.pdf", width=10, height=3)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.215, 0.145) +
  geom_hline(yintercept = 0, color="black", linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2, lineend = "round"),
        axis.ticks = element_line(linewidth = 2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()


### Extended Data Fig. 6i

####### RI2 LOS residuals 
NTTD <- df_residuals_150kb_gc[,c(5, 66)]
colnames(NTTD)[2] <- "LOS_residuals"
NTTD$treat <- "RI2"

#Violin plotting the LOS residuals for dTAG 6h
pdf("LOSr_allSPIN_global_RI2_violin_alt2.pdf", width=10, height=3)
ggplot(na.omit(NTTD), aes(x=SPIN, y=LOS_residuals, fill=SPIN)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(values = c("#9E0142",  "#F46D43", "#FDAE61", "#FFFFBF", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")) +
  ylim(-0.215, 0.145) +
  geom_hline(yintercept = 0, color="black", linetype=2, size=1.2, lineend = "round") +
  theme_classic()+
  theme(axis.text=element_text(size=15))+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  theme(axis.text=element_text(size=20)) +
  theme(axis.line = element_line(linewidth = 2, lineend = "round"),
        axis.ticks = element_line(linewidth = 2, lineend = "round"),
        axis.ticks.length = unit(0.29, "cm"))

dev.off()


###################################################  End  #####################################################################


