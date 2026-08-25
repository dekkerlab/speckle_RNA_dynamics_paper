# Liquid Chromatin Hi-C (LC-Hi-C) Analysis

## Installation

Packages needed to generate cis-percent  and LOS are included in the provided conda
environment file:

```bash
conda env create -f environment.yaml   # or: mamba env create -f environment.yaml
conda activate LOS-env
```

This installs Python, R, `r-argparse`, and (via pip) `cooler` and `bioframe`.

This environment is sufficient for the demo pipeline (cis-percent + LOS).

---

## Quick start

A small test dataset on Zenodo lets you run the core LC-Hi-C pipeline — from
balanced coolers to a genome-wide **LOS** (Loss of Structure) track — without the
full GEO deposit. The test coolers are the genome-wide K562 LC-Hi-C maps reduced
to two resolutions (50 kb, 250 kb) with their original genome-wide balancing
weights preserved, so values match the full data.

**1. Download the test data** (Zenodo DOI *10.5281/zenodo.22083200*):

```bash
mkdir -p demo/data
curl -L --fail -o demo/data/K562_notreat_LCHiC_2hPD_FA_DpnII_BR1.demo.mcool \
  "https://zenodo.org/records/22083200/files/K562_notreat_LCHiC_2hPD_FA_DpnII_BR1.demo.mcool?download=1"
curl -L --fail -o demo/data/K562_notreat_LCHiC_2hPD_FA_Mock_BR1.demo.mcool \
  "https://zenodo.org/records/22083200/files/K562_notreat_LCHiC_2hPD_FA_Mock_BR1.demo.mcool?download=1"
```

**2. cis-percent** — windowed balanced cis coverage (-r 2000000; 2Mb window) per 50 kb bin (-b 50000), 
for the Digest sample and the Mock control:

```bash
python LCHiC/cool2RangeCisPercent.py \
    -i demo/data/K562_notreat_LCHiC_2hPD_FA_DpnII_BR1.demo.mcool \
    -b 50000 -r 2000000
python LCHiC/cool2RangeCisPercent.py \
    -i demo/data/K562_notreat_LCHiC_2hPD_FA_Mock_BR1.demo.mcool \
    -b 50000 -r 2000000
```

**3. LOS** — `(Mock_cis - Digest_cis) / Mock_cis` per bin:

```bash
Rscript LCHiC/LOS.R \
    -i demo/data/K562_notreat_LCHiC_2hPD_FA_DpnII_BR1.50000_range2Mb_cispercent.bedGraph \
    -m demo/data/K562_notreat_LCHiC_2hPD_FA_Mock_BR1.50000_range2Mb_cispercent.bedGraph
```

You should get two `*_cispercent.bedGraph` files and one `*_LOS.bedGraph`
(4-column `chrom start end LOS`, higher = greater loss of local cis structure).

> **Note on short-range interactions.** No diagonals are ignored at the cooler
> level here (`cool2RangeCisPercent.py -d` defaults to 0). Instead, very
> short-range cis contacts are removed *upstream* at the pairs level during
> mapping in the `distiller-nf` pipeline, which excludes cis read pairs closer
> than 1 kb before binning. Include the following in the distiller-nf configuration
> file:
>
> ```yaml
> filters:
>     exclude1kb: '(chrom1 != chrom2) or ((chrom1 == chrom2) and (abs(pos1 - pos2) > 1e3))'
> ```

---

## Core concepts

**LOS (Loss of Structure)** is currently the main metric to assess chromtin interation 
stability using LC-Hi-C. For each bin, LOS compares the *cis* contact percentage of a 
digested (LC-Hi-C) sample to its undigested Mock control:

```
LOS = (Mock_cis - Digest_cis) / Mock_cis
```

Higher LOS = greater loss of local cis structure upon digestion (less stable
chromatin). LOS residuals are then analyzed across **SPIN states** (Speckle,
Interior, Lamina, …) and A/B compartments (PC1), with block-level aggregation to
handle spatial autocorrelation. All analysis is on **hg38**, primarily at 50 kb
and 1 Mb.

---

### The `LCHiC/` pipeline

| Step | File | Lang | Purpose |
|------|------|------|---------|
| 1 | `cool2RangeCisPercent.py` | Python | Per-bin balanced cis coverage within a ±(range/2) window, from a balanced `.cool`/`.mcool`. |
| 2 | `LOS.R` | R | LOS from a Digest cis-percent track (`-i`) and Mock control (`-m`). |
| 3 | `make_residual_table.R` | R | Joins individual LOS.bedGraph and DpnIIseq.bed files and calculates LOS residuals based on 'LOS_signal_pairs_hg38.yaml' |
| — | `LOS_signal_pairs_hg38.yaml` | config | LOS-DpnII-seq signal pairing manifest. |



---


