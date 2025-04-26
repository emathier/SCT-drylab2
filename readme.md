# Project overview
This project focuses on the single-cell multiome dataset from the study by [Azbukina, He, Lin et al.](https://doi.org/10.1101/2025.03.20.644368)
In this paper, the authors generated and analyzed paired single-cell RNA-seq (scRNA-seq) and ATAC-seq (scATAC-seq) data from human midbrain organoids at day 60. 
The goal of the study was to characterize the cellular diversity, gene regulation, and chromatin accessibility dynamics during midbrain development using an integrated single-cell multiomics approach.

# Tasks
1. Individual analysis of scRNA-seq data
  - Feature selection
  - 2D embedding (e.g., UMAP or t-SNE)
  - Cell type annotation (no integration across samples needed)

2. Individual analysis of scATAC-seq data
  - 2D embedding and visualization

3. Bi-modal integration of scRNA-seq and scATAC-seq
  - Peak-gene linkage analysis
  - Transcription factor (TF) motif enrichment analysis

4. Analysis of transcription factor regulomes and gene regulatory network (GRN) inference
  - Using the PANDO tool for GRN construction


# Organisational
The tutorial for R analysis is provided at[scMultiome_analysis_vignette.git](https://github.com/quadbio/scMultiome_analysis_vignette.git)
To get started just follow these steps:

1. Clone this repository onto your device:
   ```bash
   git clone https://github.com/emathier/SCT-drylab2.git
   cd single-cell-multiome-organoids

2. Install the necessary R packages as described in the tutorial.
3. Dont forget to save (= commit and push the changes) to this repository with
   ```bash
    git add .
    git commit -m "Add initial README with project overview and goals"
    git push origin main
    git status
