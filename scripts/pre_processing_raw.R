library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
library(UpSetR)
theme_set(theme_cowplot())

# Directory names for all filtered data
data_dirs <- c(
  "~/Research_Project/data/raw/2018_control_leukocytes", 
  "~/Research_Project/data/raw/2018_control_mesenchyme",
  "~/Research_Project/data/raw/2018_T24_repairing_leukocytes", 
  "~/Research_Project/data/raw/2018_T24_repairing_mesenchyme",
  "~/Research_Project/data/raw/2020_control_leukocytes", 
  "~/Research_Project/data/raw/2020_control_mesenchyme",
  "~/Research_Project/data/raw/2020_T48_repairing_leukocytes", 
  "~/Research_Project/data/raw/2020_T48_repairing_mesenchyme"
)

# Named list to store Seurat objects
seurat_objects <- list()

# Loop through each directory and create Seurat objects
for (dir in data_dirs) {
  # Extract a short label from the folder name
  label <- basename(dir)
  
  # Read using ReadMtx rather than 10x
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
  seurat_obj$timepoint <- parts[2]
  seurat_obj$celltype <- parts[3]
  
  
  # Store in list
  seurat_objects[[label]] <- seurat_obj
}

# Merge controls by cell type across years
controls_leukocytes <- merge(
  seurat_objects[["2018_control_leukocytes"]],
  y = seurat_objects[["2020_control_leukocytes"]],
  add.cell.ids = c("2018_leuko", "2020_leuko"),
  project = "controls_leukocytes"
)
controls_mesenchyme <- merge(
  seurat_objects[["2018_control_mesenchyme"]],
  y = seurat_objects[["2020_control_mesenchyme"]],
  add.cell.ids = c("2018_mesen", "2020_mesen"),
  project = "controls_mesenchyme"
)

# Store merged controls in list, drop originals
seurat_objects[["controls_leukocytes"]] <- controls_leukocytes
seurat_objects[["controls_mesenchyme"]] <- controls_mesenchyme

seurat_objects[["2018_control_leukocytes"]]   <- NULL
seurat_objects[["2020_control_leukocytes"]]   <- NULL
seurat_objects[["2018_control_mesenchyme"]]   <- NULL
seurat_objects[["2020_control_mesenchyme"]]   <- NULL

# Combine all 6 Seurat objects into a big object
big_obj <- merge(x = seurat_objects[["controls_leukocytes"]],
                 y = list(
                   seurat_objects[["controls_mesenchyme"]],
                   seurat_objects[["2018_T24_repairing_leukocytes"]], 
                   seurat_objects[["2018_T24_repairing_mesenchyme"]],
                   seurat_objects[["2020_T48_repairing_leukocytes"]], 
                   seurat_objects[["2020_T48_repairing_mesenchyme"]]
                 ),
                 add.cell.ids = c("cont_leuk", "cont_mes", "24_leuk", "24_mes", "48_leuk", "48_mes"),
                 project = "all_expmnts")

big_obj$short_id <- recode(big_obj$orig.ident,
                           "2018_control_leukocytes" = "18_cont_leuk",
                           "2018_control_mesenchyme" = "18_cont_mes",
                           "2018_T24_repairing_leukocytes" = "24_leuk",
                           "2018_T24_repairing_mesenchyme" = "24_mes",
                           "2020_T48_repairing_leukocytes" = "48_leuk",
                           "2020_T48_repairing_mesenchyme" = "48_mes",
                           "2020_control_leukocytes" = "20_cont_leuk",
                           "2020_control_mesenchyme" = "20_cont_mes")

# Calculate percent.mt and percent.ribo
big_obj[["percent.mt"]] <- PercentageFeatureSet(big_obj, pattern = "^mt-")
big_obj[["percent.ribo"]] <- PercentageFeatureSet(big_obj, pattern = "^Rpl|^Rps")

qc_summary <- data.frame(
  total_cells = ncol(big_obj),
  cells_low_features = sum(big_obj$nFeature_RNA < 200),
  cells_high_features = sum(big_obj$nFeature_RNA > 5000),
  cells_high_mito = sum(big_obj$percent.mt >= 10),
  cells_high_ribo = sum(big_obj$percent.ribo >= 10)
)
qc_summary
