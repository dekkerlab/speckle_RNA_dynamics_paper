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

Typical installation of software and dependencies should take a few minutes. Run time for figure generation 
depends on the analysis. While most plots can be generated within a few minutes, operations involving scaling
or calculating expected interactions can take 10-15minutes or more. All figures generated in **R Studio** list 
their required packages at the top of each script, so no separate R environment file is provided — install the 
packages named in the relevant `Fig_*/*.R` script.

---

As shown in the **LCHiC** directory, a small test dataset on Zenodo lets you run the core LC-Hi-C pipeline to 
generate LOS and LOS residuals from balanced coolers. Most heatmaps and stacks can be plotted directly 
using these test coolers. 

### Figure sub-directories (brief)

- **Fig_1** — compartments/PC1, raw LOS tracks (incl. heat-shock), Hi-C heatmaps.
- **Fig_2** — Crane maps; SON immunofluorescence intensity from `.lif` imaging (Aivia output consolidated in R).
- **Fig_3** — DNA-FISH: nucleus/foci segmentation (`scikit-image`) and SON-proximity quantification.
- **Fig_4** — SLAM-seq nascent transcription (all/exons/introns) vs. structure.
- **Fig_5** — differential expression integrated with the 1 Mb residual table.
- **Fig_6** — cCRE and scaled gene-body pileups/stackups over LOS bigWigs (`bbi.stackup`).
- **ED_Fig_5** — GO enrichment (`clusterProfiler`).

> **Paths and File names** many scripts contain hard-coded absolute paths and alternate file names; edit them to point to your
> own data before rerunning.
