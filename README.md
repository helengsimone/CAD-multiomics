# Multi-omics analysis reveals a non-coding variant connecting CRHR1 signalling to coronary artery disease

This repository contains the custom R code used for GWAS variant prioritisation via sequential epigenomic integration, as described in the manuscript.

---

## Script

`post_gwas_pipeline.R`

---

## Dependencies

R packages:

- `data.table`
- `rtracklayer` (Bioconductor)
- `GenomicRanges` (Bioconductor)

Install with:

```r
install.packages("data.table")
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("rtracklayer", "GenomicRanges"))
```

---

## Input files

| File | Source | Description |
|------|--------|-------------|
| `GWAS.csv` | CARDIoGRAM consortium | GWAS summary statistics (hg18) |
| `hg18ToHg19.over.chain.gz` | UCSC Genome Browser | LiftOver chain file |
| `dna_seq.bed` | ENCODE / GEO | DNase-seq peaks, human aortic adventitial fibroblasts |
| `GSM906392_UCSD.Aorta.H3K27ac.STL003.bed` | GEO: GSM906392 | H3K27ac ChIP-seq fragments, human aorta |
| `GSM910566_UCSD.Aorta.H3K36me3.STL003.bed` | GEO: GSM910566 | H3K36me3 ChIP-seq fragments, human aorta |

---

## Pipeline

The script runs five sequential filtering steps:

1. **LiftOver** — GWAS coordinates converted from hg18 to hg19
2. **Nominal p-value filter** — variants with p < 0.05 retained
3. **DNase-seq overlap** — variants intersected with open chromatin regions
4. **H3K27ac overlap** — variants intersected with active regulatory elements
5. **H3K36me3 overlap** — variants intersected with actively transcribed regions

---

## Output

`CAD_p005_DNase_H3K27ac_H3K36me3.csv` — final prioritised variant set passing all epigenomic filters

---

## Reference genome

All analyses use **GRCh37/hg19**.
