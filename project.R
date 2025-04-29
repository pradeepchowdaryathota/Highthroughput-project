# 1. Load Required Libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)


setwd("/Users/pradeepchowdary/Desktop/Hightroughput/GSE180286_RAW")

# list all the expression matrix files
sample_files <- list.files(pattern = "expression_matrix\\.txt$")

# Read each file and create a Seurat object
seurat_list <- lapply(sample_files, function(file) {
  counts <- read.table(file, header = TRUE, row.names = 1, sep = "\t")
  CreateSeuratObject(counts = counts, project = gsub("_expression_matrix.txt", "", file))
})

# Merge all Seurat objects
seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = gsub("_expression_matrix.txt", "", sample_files))

# Save merged object
saveRDS(seurat_obj, "breast_ln_merged_seurat.rds")


# Quality Control
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

# Normalize and Scale Data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# Dimensionality Reduction and Clustering
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Plot UMAP
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("UMAP Plot: Breast Cancer Lymph Node Metastasis Cells")

# Define Checkpoint Genes
checkpoint_genes <- c("PDCD1", "CTLA4", "TIGIT", "SIGLEC10", "CD274", "LAG3", "HAVCR2", "CD52")
present_genes <- intersect(checkpoint_genes, rownames(seurat_obj))

# Dot Plot
DotPlot(seurat_obj, features = present_genes, group.by = "seurat_clusters") +
  RotatedAxis() + ggtitle("Checkpoint Gene Expression Across Clusters")

# Violin Plot
VlnPlot(
  seurat_obj,
  features = present_genes,
  group.by = "seurat_clusters",
  pt.size = 0.1,           
  ncol = 2                 
) +
  ggtitle("Improved Violin Plot: Checkpoint Gene Expression per Cluster") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_brewer(palette = "Paired")  

# Export Average Expression Table
avg_expr <- AverageExpression(seurat_obj, features = present_genes, group.by = "seurat_clusters", return.seurat = FALSE)
write.csv(avg_expr$RNA, "Table1_MeanExpression_Checkpoints.csv")
