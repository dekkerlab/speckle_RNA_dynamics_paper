# library(tidyverse)
library(stringr)
library(yaml)


resolution <- "50kb"

data_folder <- paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/",resolution,"/", sep="")





bins <- read.table(paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Sequencing/data/binning/hg38_", resolution,".bed",sep=""),
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end"))
bins <- bins[ - grep("chrY", bins$chrom),]
bins <- bins[ - grep("chrM", bins$chrom),]

bins$chrom <- factor(bins$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
bins <- bins[order(bins$chrom, bins$start,bins$end),]
rownames(bins) <- NULL


eigen <- read.table("C:/Users/lafon/UMass Medical School Dropbox/Denis Lafontaine/UmassMed/Dekker/Sequencing/Analysis/20250129_Mega_maps_hg38/analysis/eigen/K562_hg38.50kb.eigs.cis.vecs.txt",
                   sep="\t", header=TRUE)
eigen <- eigen[ - grep("chrY", eigen$chrom),]
eigen <- eigen[ - grep("chrM", eigen$chrom),]
colnames(eigen)[4] <- "PC1"
eigen[is.na(eigen$PC1),"PC1"] <- NA


SPIN <- read.table(paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Sequencing/data/K562_data/SPIN/hg38/SPIN_", resolution,"_rebinned.bed",sep=""), sep="\t",
                   header=TRUE)
SPIN <- SPIN[, !names(SPIN) %in% c("n_bins", "n_unique")]
SPIN$SPIN <- factor(SPIN$SPIN, order = TRUE, levels =c("Speckle", "Interior_Act1", "Interior_Act2", "Interior_Act3", "Interior_Repr1","Interior_Repr2","Near_Lm1","Near_Lm2","Lamina","Lamina_Like"))



############### Functions  #################################
### Extract column name from sample name


col_name_LOS <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  LOS_name <- paste("LOS_",unlist(str_split(b,"_"))[4], "_", unlist(str_split(b,"-"))[4],"_",unlist(str_split(b,"-"))[5],"_",unlist(str_split(b,"-"))[6],"_", unlist(str_split(b,"-"))[7],"_",unlist(str_split(b,"-"))[2], sep="")
}


col_name_DpnIIseq <- function(x){
  a <- unlist(str_split(x,"/"))  
  b <- tail(a,1)
  DpnII_name <- paste("signal_", unlist(str_split(b,"-"))[4],"_",unlist(str_split(b,"-"))[5],"_",unlist(str_split(b,"-"))[6],"_",unlist(str_split(b,"-"))[7],"_",unlist(str_split(b,"-"))[2],sep = "")
}

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

#### list data from each folder (cis_coverage, cis_percent, LOS, smooth_LOS)
#cis_percentList <- list.files(paste(data_folder,"cis_percent",sep=""), pattern = "cispercent.bedGraph", full.name = TRUE)

#LOS_coverageList <- list.files(paste(data_folder,"cis_coverage",sep=""), pattern = "cov_LOS", full.name = TRUE)

LOSList <- list.files(paste(data_folder, "/LOS", sep =""), pattern = "LOS", full.name = TRUE)

#smooth_LOSList <- list.files(paste(data_folder,"LOSsm", sep=""), pattern = "smooth.bedGraph", full.name = TRUE)

DpnIIList <- list.files(paste(data_folder,"DpnII_seq", sep=""), pattern = "copy_correct_coverage_", full.name = TRUE)



LOS_names <- lapply(LOSList, col_name_LOS)
LOS_names <- unlist(LOS_names) 
LOS_Data <- lapply(1:length(LOSList), function(x) { read.delim(LOSList[[x]], header = TRUE)})


for (i in c(1:length(LOS_Data))) {
  names(LOS_Data[[i]]) <- c("chrom","start","end",LOS_names[[i]])
}
names(LOS_Data) <- c(LOS_names)

for (i in c(1:length(LOS_Data))) {
  LOS_Data[[i]][2] <- LOS_Data[[i]][2] + 1
}

### Remove chromosomes Y and M from LOS data

for (i in c(1:length(LOS_Data))) {
  LOS_Data[[i]][!grepl("chrY", LOS_Data[[i]]$chrom),]
}



#### Load and munge DpnII-seq tables
signal_names <- lapply(DpnIIList, col_name_DpnIIseq)
signal_names <- unlist(signal_names)
signal_Data <- lapply(1:length(DpnIIList), function(x) { read.delim(DpnIIList[[x]], header = FALSE,col.names= c("chrom", "start", "end", "signal"))})
for (i in c(1:length(signal_Data))) {
  names(signal_Data[[i]]) <- c("chrom","start","end",signal_names[[i]])
}
names(signal_Data) <- c(signal_names)


#### Make d table ####

d <- data.frame(eigen[,1:4],SPIN["SPIN"])


for (i in c(1:length(signal_Data))) {
  d <- merge(d, signal_Data[[i]], by=c("chrom", "start", "end"), all.x= TRUE)
}


for (i in c(1:length(LOS_Data))) {
  d <- merge(d, LOS_Data[[i]], by=c("chrom", "start", "end"), all.x= TRUE)
}


d$chrom <- factor(d$chrom, levels=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))
d <- d[order(d$chrom, d$start,d$end),]
rownames(d) <- NULL



setwd(paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/",resolution,"/dataframes/", sep=""))
write.table(d, paste("d_20260812_",resolution,".bed",sep=""), col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)




df_residuals <- d

# Read the YAML configuration file
config <- yaml.load_file(paste0("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/config/", resolution,"/LOS_signal_pairs_hg38.yaml"))
setwd(paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/",resolution,"/LOS_residuals",sep=""))
# Loop through the key-value pairs in the YAML file
for (pair in names(config)) {
  LOS_name <- config[[pair]]$LOS
  signal_name <- config[[pair]]$signal
  #LOS_name <- config[[1]]$LOS
  #signal_name <- config[[1]]$signal  
  
  d_clean <-  na.omit(df_residuals[, c("chrom", "start", "end", LOS_name, signal_name)])
  
  los_sample <- d_clean[,LOS_name]
  signal_sample <- d_clean[,signal_name]

  n <- mean(na.omit(d_clean[,signal_name]))/4
  
  ma <- get_moving_average(signal_sample, los_sample, n, 1, 0)
  los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", signal_name, LOS_name)], ma, 0)
  names(los_r)[4] <- paste("LOS_residuals_", unlist(str_split(LOS_name,"_"))[2], "_", unlist(str_split(LOS_name,"_"))[3], "_",unlist(str_split(LOS_name,"_"))[4], "_",unlist(str_split(LOS_name,"_"))[5], "_",unlist(str_split(LOS_name,"_"))[6], "_",tail(unlist(str_split(LOS_name,"_"))[-1], n=1), sep="")
  
  df_residuals <- merge(df_residuals,los_r, by=c("chrom", "start", "end"), all.x= TRUE)
  los_r_clean <- merge(bins, los_r, by=c("chrom", "start", "end"), all.x= TRUE)
  write.table(los_r_clean, paste("LOS_residuals_", unlist(str_split(LOS_name,"_"))[2], "_", unlist(str_split(LOS_name,"_"))[3], "_",unlist(str_split(LOS_name,"_"))[4], "_",unlist(str_split(LOS_name,"_"))[5], "_",unlist(str_split(LOS_name,"_"))[6],"_",tail(unlist(str_split(LOS_name,"_"))[-1], n=1),".bedGraph", sep=""), col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)
}

setwd(paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/",resolution,"/dataframes/", sep=""))
write.table(df_residuals, paste("df_residuals_20260604_",resolution,".bed",sep=""), col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)




#######  Add diff columns  ####

####Calculate differences in LOS residuals between appropriate samples

df_residuals["LOS_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif"] <- df_residuals["LOS_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103"] - df_residuals["LOS_range2Mb_K562_PIIBTD240mPD_DpnII_R1_20200103"]
df_residuals["LOS_range2Mb_K562_PIIBNTTD_Mock_R2_20201102_dif"] <- df_residuals["LOS_range2Mb_K562_PIIBNT_Mock_R2_20201102"] - df_residuals["LOS_range2Mb_K562_PIIBTD_DpnII_R2_20201102"]
df_residuals["LOS_range2Mb_HAP1RPB1AIDGFP_IAA240mPD_DpnII_R1_20210721_dif"] <- df_residuals["LOS_range2Mb_HAP1RPB1AIDGFP_WAT240mPD_DpnII_R1_20210721"] - df_residuals["LOS_range2Mb_HAP1RPB1AIDGFP_IAA240mPD_DpnII_R1_20210721"]
df_residuals["LOS_range2Mb_K562_RNaseA4hPD_DpnII_R1_20220120_dif"] <- df_residuals["LOS_range2Mb_K562_RNaseIN4hPD_DpnII_R1_20220120"] - df_residuals["LOS_range2Mb_K562_RNaseA4hPD_DpnII_R1_20220120"]
df_residuals["LOS_range2Mb_K562_U1AMO4hPD_DpnII_R1_20220922_dif"] <- df_residuals["LOS_range2Mb_K562_contAMO4hPD_DpnII_R1_20220922"] - df_residuals["LOS_range2Mb_K562_U1AMO4hPD_DpnII_R1_20220922"]
df_residuals["LOS_range2Mb_K562_PIIBDRB_DpnII_R1_20230223_dif"] <- df_residuals["LOS_range2Mb_K562_PIIBNT_DpnII_R3_20230223"] - df_residuals["LOS_range2Mb_K562_PIIBDRB_DpnII_R1_20230223"]
df_residuals["LOS_range2Mb_K562_PIIBTPL_DpnII_R1_20230223_dif"] <- df_residuals["LOS_range2Mb_K562_PIIBNT_DpnII_R3_20230223"] - df_residuals["LOS_range2Mb_K562_PIIBTPL_DpnII_R1_20230223"]
df_residuals["LOS_range2Mb_K562_PIIBNTTD_DpnII_R3_20230223_dif"] <- df_residuals["LOS_range2Mb_K562_PIIBNT_DpnII_R3_20230223"] - df_residuals["LOS_range2Mb_K562_PIIBTD_DpnII_R3_20230223"]
df_residuals["LOS_range2Mb_K562_HS42240mPD_DpnII_R1_20230223_dif"] <- df_residuals["LOS_range2Mb_K562_HS37240mPD_DpnII_R1_20230223"] - df_residuals["LOS_range2Mb_K562_HS42240mPD_DpnII_R1_20230223"]
df_residuals["LOS_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628_dif"] <- df_residuals["LOS_range2Mb_K562_RNaseIN4hPD_DpnII_R2_20230628"] - df_residuals["LOS_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628"]
df_residuals["LOS_range2Mb_K562_HS42120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_range2Mb_K562_HS37120mPD_DpnII_R2_20231218"] - df_residuals["LOS_range2Mb_K562_HS42120mPD_DpnII_R2_20231218"]
df_residuals["LOS_range2Mb_K562_old120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_range2Mb_K562_new120mPD_DpnII_R2_20231218"] - df_residuals["LOS_range2Mb_K562_old120mPD_DpnII_R2_20231218"]
df_residuals["LOS_range2Mb_K562_U1AMO120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_range2Mb_K562_contAMO120mPD_DpnII_R2_20231218"] - df_residuals["LOS_range2Mb_K562_U1AMO120mPD_DpnII_R2_20231218"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_6hdTAGDMSO120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_DMSO120mPD_DpnII_R2_20231218"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAG120mPD_DpnII_R2_20231218"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_Neg120mPD_DpnII_R2_20231218"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAG120mPD_DpnII_R2_20231218"]
df_residuals["LOS_range2Mb_K562_CdSulfate_120mPD_DpnII_20240327_dif"] <- df_residuals["LOS_range2Mb_K562_contCD_120mPD_DpnII_20240327"] - df_residuals["LOS_range2Mb_K562_CdSulfate_120mPD_DpnII_20240327"]
df_residuals["LOS_range2Mb_K562_old120mPD_DpnII_R2_20240327_dif"] <- df_residuals["LOS_range2Mb_K562_contCD_120mPD_DpnII_20240327"] - df_residuals["LOS_range2Mb_K562_old120mPD_DpnII_R2_20240327"]
df_residuals["LOS_range2Mb_K562_RNaseA120hPD_DpnII_R3_20240611_dif"] <- df_residuals["LOS_range2Mb_K562_RNaseIN120hPD_DpnII_R3_20240611"] - df_residuals["LOS_range2Mb_K562_RNaseA120hPD_DpnII_R3_20240611"]
df_residuals["LOS_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif"] <- df_residuals["LOS_range2Mb_K562_contAMO120mPD_DpnII_R4_20240828"] - df_residuals["LOS_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_11hdTAG120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGDMSO120mPD_DpnII_R1_20240828"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_8hdTAG120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_NegHS37120mPD_DpnII_R1_20240828"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGHS37120mPD_DpnII_R1_20240828"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_NegNTTD120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_NegTD120mPD_DpnII_R1_20240828"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGNTTD120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGDMSO120mPD_DpnII_R1_20240828"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGTD120mPD_DpnII_R1_20240828"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_NegHS42120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_NegHS37120mPD_DpnII_R1_20240828"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_NegHS42120mPD_DpnII_R1_20240828"]
df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGHS42120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGHS37120mPD_DpnII_R1_20240828"] - df_residuals["LOS_range2Mb_K562SpeckledTAG_dTAGHS42120mPD_DpnII_R1_20240828"]
df_residuals["LOS_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif"] <- df_residuals["LOS_range2Mb_K562HS37_120mPD_DpnII_T3_20250624"] - df_residuals["LOS_range2Mb_K562HS42_120mPD_DpnII_T3_20250624"]

df_residuals["LOS_residuals_range2Mb_K562_PIIBNTTD240mPD_DpnII_R1_20200103_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_PIIBNT240mPD_DpnII_R1_20200103"] - df_residuals["LOS_residuals_range2Mb_K562_PIIBTD240mPD_DpnII_R1_20200103"]
df_residuals["LOS_residuals_range2Mb_K562_PIIBNTTD_Mock_R2_20201102_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_PIIBNT_Mock_R2_20201102"] - df_residuals["LOS_residuals_range2Mb_K562_PIIBTD_DpnII_R2_20201102"]
df_residuals["LOS_residuals_range2Mb_HAP1RPB1AIDGFP_IAA240mPD_DpnII_R1_20210721_dif"] <- df_residuals["LOS_residuals_range2Mb_HAP1RPB1AIDGFP_WAT240mPD_DpnII_R1_20210721"] - df_residuals["LOS_residuals_range2Mb_HAP1RPB1AIDGFP_IAA240mPD_DpnII_R1_20210721"]
df_residuals["LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R1_20220120_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_RNaseIN4hPD_DpnII_R1_20220120"] - df_residuals["LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R1_20220120"]
df_residuals["LOS_residuals_range2Mb_K562_U1AMO4hPD_DpnII_R1_20220922_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_contAMO4hPD_DpnII_R1_20220922"] - df_residuals["LOS_residuals_range2Mb_K562_U1AMO4hPD_DpnII_R1_20220922"]
df_residuals["LOS_residuals_range2Mb_K562_PIIBDRB_DpnII_R1_20230223_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_PIIBNT_DpnII_R3_20230223"] - df_residuals["LOS_residuals_range2Mb_K562_PIIBDRB_DpnII_R1_20230223"]
df_residuals["LOS_residuals_range2Mb_K562_PIIBTPL_DpnII_R1_20230223_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_PIIBNT_DpnII_R3_20230223"] - df_residuals["LOS_residuals_range2Mb_K562_PIIBTPL_DpnII_R1_20230223"]
df_residuals["LOS_residuals_range2Mb_K562_PIIBNTTD_DpnII_R3_20230223_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_PIIBNT_DpnII_R3_20230223"] - df_residuals["LOS_residuals_range2Mb_K562_PIIBTD_DpnII_R3_20230223"]
df_residuals["LOS_residuals_range2Mb_K562_HS42240mPD_DpnII_R1_20230223_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_HS37240mPD_DpnII_R1_20230223"] - df_residuals["LOS_residuals_range2Mb_K562_HS42240mPD_DpnII_R1_20230223"]
df_residuals["LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_RNaseIN4hPD_DpnII_R2_20230628"] - df_residuals["LOS_residuals_range2Mb_K562_RNaseA4hPD_DpnII_R2_20230628"]
df_residuals["LOS_residuals_range2Mb_K562_HS42120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_HS37120mPD_DpnII_R2_20231218"] - df_residuals["LOS_residuals_range2Mb_K562_HS42120mPD_DpnII_R2_20231218"]
df_residuals["LOS_residuals_range2Mb_K562_old120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_new120mPD_DpnII_R2_20231218"] - df_residuals["LOS_residuals_range2Mb_K562_old120mPD_DpnII_R2_20231218"]
df_residuals["LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_contAMO120mPD_DpnII_R2_20231218"] - df_residuals["LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R2_20231218"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_6hdTAGDMSO120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_DMSO120mPD_DpnII_R2_20231218"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAG120mPD_DpnII_R2_20231218"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_6hdTAGNeg120mPD_DpnII_R2_20231218_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_Neg120mPD_DpnII_R2_20231218"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAG120mPD_DpnII_R2_20231218"]
df_residuals["LOS_residuals_range2Mb_K562_CdSulfate_120mPD_DpnII_20240327_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_contCD_120mPD_DpnII_20240327"] - df_residuals["LOS_residuals_range2Mb_K562_CdSulfate_120mPD_DpnII_20240327"]
df_residuals["LOS_residuals_range2Mb_K562_old120mPD_DpnII_R2_20240327_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_contCD_120mPD_DpnII_20240327"] - df_residuals["LOS_residuals_range2Mb_K562_old120mPD_DpnII_R2_20240327"]
df_residuals["LOS_residuals_range2Mb_K562_RNaseA120hPD_DpnII_R3_20240611_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_RNaseIN120hPD_DpnII_R3_20240611"] - df_residuals["LOS_residuals_range2Mb_K562_RNaseA120hPD_DpnII_R3_20240611"]
df_residuals["LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562_contAMO120mPD_DpnII_R4_20240828"] - df_residuals["LOS_residuals_range2Mb_K562_U1AMO120mPD_DpnII_R4_20240828"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_11hdTAG120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGDMSO120mPD_DpnII_R1_20240828"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_8hdTAG120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegHS37120mPD_DpnII_R1_20240828"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGHS37120mPD_DpnII_R1_20240828"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegNTTD120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R1_20240828"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegTD120mPD_DpnII_R1_20240828"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGNTTD120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGDMSO120mPD_DpnII_R1_20240828"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGTD120mPD_DpnII_R1_20240828"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegHS42120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegHS37120mPD_DpnII_R1_20240828"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegHS42120mPD_DpnII_R1_20240828"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGHS42120mPD_DpnII_R1_20240828_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGHS37120mPD_DpnII_R1_20240828"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGHS42120mPD_DpnII_R1_20240828"]
df_residuals["LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624_dif"] <- df_residuals["LOS_residuals_range2Mb_K562HS37_120mPD_DpnII_T3_20250624"] - df_residuals["LOS_residuals_range2Mb_K562HS42_120mPD_DpnII_T3_20250624"]

df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_6hdTAG120mPD_DpnII_R1_20260616_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_Neg120mPD_DpnII_R2_20260616"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAG120mPD_DpnII_R2_20260616"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_10hdTAG120mPD_DpnII_R1_20260616_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R2_20260616"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGDMSO120mPD_DpnII_R2_20260616"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegNTTD120mPD_DpnII_R1_20260616_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegDMSO120mPD_DpnII_R2_20260616"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_NegTD120mPD_DpnII_R2_20260616"]
df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGNTTD120mPD_DpnII_R1_20260616_dif"] <- df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGDMSO120mPD_DpnII_R2_20260616"] - df_residuals["LOS_residuals_range2Mb_K562SpeckledTAG_dTAGTD120mPD_DpnII_R2_20260616"]



setwd(paste("C:/Users/lafon/Dropbox (UMass Medical School)/UmassMed/Dekker/Manuscripts/2022_RNA_Genome Stability/analysis/data/",resolution,"/dataframes/", sep=""))
write.table(df_residuals, paste("df_residuals_20260812_",resolution,".bed",sep=""), col.names = TRUE, row.names = FALSE, sep = '\t',quote = FALSE)

