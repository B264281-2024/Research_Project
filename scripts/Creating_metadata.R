#Initial data analysis - BC Research Project
#date: 26.5.2025

#setup required packages
#BiocManager::install("Seurat")
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
library(UpSetR)
theme_set(theme_cowplot())

#Directory names for all filtered data
data_dirs <- c(
  "~/Research_Project/data/filtered/2018_control_leukocytes", 
  "~/Research_Project/data/filtered/2018_control_mesenchyme",
  "~/Research_Project/data/filtered/2018_T24_repairing_leukocytes", 
  "~/Research_Project/data/filtered/2018_T24_repairing_mesenchyme",
  "~/Research_Project/data/filtered/2020_control_leukocytes", 
  "~/Research_Project/data/filtered/2020_control_mesenchyme",
  "~/Research_Project/data/filtered/2020_T48_repairing_leukocytes", 
  "~/Research_Project/data/filtered/2020_T48_repairing_mesenchyme"
)

#setwd("~/Research_Project/data/filtered/2018_control_leukocytes")
#readLines("barcodes.tsv", n = 5)
#read.table("features.tsv", sep = "\t", head = FALSE, nrows = 5)

# Named list to store Seurat objects
seurat_objects <- list()

# Loop through each directory and create Seurat objects
for (dir in data_dirs) {
  # Extract a short label from the folder name
  label <- basename(dir)
  
  # Read using ReadMtx
  data <- ReadMtx(
    mtx = file.path(dir, "matrix.mtx.gz"),
    features = file.path(dir, "features.tsv.gz"),
    cells = file.path(dir, "barcodes.tsv.gz"),
    feature.column = 2
  )
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = data, project = label)
  parts <- strsplit(label, "_")[[1]]
  seurat_obj$year <- parts[1]
  seurat_obj$condition <- parts[2]
  seurat_obj$tissue <- parts[3]
  
  
  # Store in list
  seurat_objects[[label]] <- seurat_obj
}

# Merge leukocyte controls
#merged_leukocytes <- merge(
#  seurat_objects[["2018_control_leukocytes"]],
#  y = seurat_objects[["2020_control_leukocytes"]],
#  add.cell.ids = c("2018_leuko", "2020_leuko"),
#  project = "Control_Leukocytes"
#)

# Merge mesenchyme controls
#merged_mesenchyme <- merge(
#  seurat_objects[["2018_control_mesenchyme"]],
#  y = seurat_objects[["2020_control_mesenchyme"]],
#  add.cell.ids = c("2018_mes", "2020_mes"),
#  project = "Control_Mesenchyme"
#)

# Create new list of processed objects
#processed_objects <- list()
#processed_objects[["control_leukocytes"]] <- merged_leukocytes
#processed_objects[["control_mesenchyme"]] <- merged_mesenchyme

# Add non-control datasets
#non_control_labels <- setdiff(names(seurat_objects), c(
#  "2018_control_leukocytes",
#  "2020_control_leukocytes",
#  "2018_control_mesenchyme",
#  "2020_control_mesenchyme"
#))
#processed_objects <- c(processed_objects, seurat_objects[non_control_labels])


for (label in names(seurat_objects)) {
  obj <- seurat_objects[[label]]
  
  # Add optional QC metrics
  #obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  #obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")
  
  # Standard Seurat preprocessing
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 20)
  obj <- RunUMAP(obj, reduction="pca", dims=1:20)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.2)
  
  # Store back into the list
  seurat_objects[[label]] <- obj
}

#VlnPlot(object, features = c("nFeature_RNA", "percent.mito"), group.by = "seurat_clusters")

# Generate UMAPs
umap_plots <- lapply(names(seurat_objects), function(label) {
  DimPlot(seurat_objects[[label]], reduction = "umap", label = TRUE) + 
    ggtitle(label)
})

# Display all UMAPs
batch_size <- 4
n <- length(umap_plots)
for (i in seq(1, n, by = batch_size)) {
  batch <- umap_plots[i:min(i + batch_size - 1, n)]
  print(plot_grid(plotlist = batch, ncol = 2))
}


