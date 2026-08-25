# Figures directory

Two additional conda environments were constructed **for reproducing figures only**:

1. **Genomic analysis & plotting** (environment_open2C.yaml) — for the Jupyter notebooks that analyze and
   plot all genomic data (heatmaps, Crane/saddle maps, cCRE and gene-body
   pileups).
2. **Imaging analysis** (environment_imaging.yaml) — for the DNA-FISH / microscopy pipelines (`.lif`
   reading, segmentation, intensity quantification).

## Installation

```bash
conda env create -f <enter appropriate environment.yaml>   
conda activate <env_name>
```

All figures generated in **R** list their required packages at the top of each script, so no separate R 
environment file is provided — install the packages named in the relevant `Fig_*/*.R` script.

---

As shown in the **LCHiC** directory, a small test dataset on Zenodo lets you run the core LC-Hi-C pipeline — from
balanced coolers to a genome-wide **LOS** (Loss of Structure) track — without the full GEO deposit. The test 
coolers are the genome-wide K562 LC-Hi-C maps reduced to two resolutions (50 kb, 250 kb) with their original 
genome-wide balancing weights preserved, so values match the full data. Heatmaps and stacks can be plotted directly 
using these test coolers. 


**1. Download the test data** (Zenodo DOI *(add concept DOI)*):

```bash
mkdir -p demo/data
curl -L --fail -o demo/data/K562_notreat_LCHiC_2hPD_FA_DpnII_BR1.demo.mcool \
  "https://zenodo.org/records/<RECORD_ID>/files/K562_notreat_LCHiC_2hPD_FA_DpnII_BR1.demo.mcool?download=1"
curl -L --fail -o demo/data/K562_notreat_LCHiC_2hPD_FA_Mock_BR1.demo.mcool \
  "https://zenodo.org/records/<RECORD_ID>/files/K562_notreat_LCHiC_2hPD_FA_Mock_BR1.demo.mcool?download=1"
```



### Figure sub-directories (brief)

- **Fig_1** — compartments/PC1, raw LOS tracks (incl. heat-shock), Hi-C heatmaps.
- **Fig_2** — Crane maps; SON immunofluorescence intensity from `.lif` imaging (Aivia output consolidated in R).
- **Fig_3** — DNA-FISH: nucleus/foci segmentation (`scikit-image`) and SON-proximity quantification.
- **Fig_4** — SLAM-seq nascent transcription (all/exons/introns) vs. structure.
- **Fig_5** — differential expression integrated with the 1 Mb residual table.
- **Fig_6** — cCRE and scaled gene-body pileups/stackups over LOS bigWigs (`bbi.stackup`).
- **ED_Fig_5** — GO enrichment (`clusterProfiler`).

> **Paths:** many scripts contain hard-coded absolute paths to the original
> analysis environment; edit them to point at your own data before rerunning.
