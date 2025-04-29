# Highthroughput-project

# Breast Cancer scRNA-seq Checkpoint Gene Analysis

This project analyzes single-cell RNA sequencing (scRNA-seq) data from breast cancer lymph node metastases using publicly available data from GEO (GSE180286). The focus is on exploring immune checkpoint gene expression across different cell clusters.

---

## Overview

- **Dataset**: GSE180286 (breast cancer primary + lymph node scRNA-seq)
- **Objective**: Identify immune checkpoint expression patterns in metastatic lymph node cell clusters
- **Tools used**: R, Seurat, ggplot2

---

## Key Steps

1. **Data Download & Preprocessing**  
   Extracted count matrices from GSE180286 and merged into a single Seurat object.

2. **Clustering & Visualization**  
   Performed normalization, PCA, clustering, and UMAP visualization.

3. **Checkpoint Gene Analysis**  
   Focused on genes like `PDCD1`, `CTLA4`, `TIGIT`, `CD274`, etc.  
   Created UMAP, DotPlot, and Violin Plot visualizations.

4. **Exported Summary Table**  
   Average expression values for checkpoint genes across clusters.

---

## Output Files

- `Figure1_UMAP_CellTypes.png` – UMAP showing clustered cells  
- `Figure2_Dotplot_Checkpoints.png` – Dot plot of checkpoint genes  
- `Figure3_Violin_Checkpoints.png` – Violin plot showing expression distribution  
- `Table1_MeanExpression_Checkpoints.csv` – Expression table of selected genes

---

## How to Run

1. Download GSE180286 and extract into your working folder  
2. Open R and run the provided script:  
   ```r
   source("scripts/analysis.R")
