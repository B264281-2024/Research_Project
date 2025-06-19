#Analysis of single cell RNA seq data using Seurat
#Naive and repairing leukocyte datasets-UMAP analysis updated to include data from 2018&2020
#date: 10/08/2021

#setup required packages
BiocManager::install("Seurat")
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
library(UpSetR)
theme_set(theme_cowplot())
source("SCRIPTS/Utils.R")

setwd("/media/pkirkwoo/Data/SCSeq_Analysis/Analysis_1")

#=================================================================================================================================================================================================================
#read in data to create separate seurat objects 
data <- Read10X(data.dir = "RAW_FILES/2018_T24_repairing_leukocytes/outs/filtered_feature_bc_matrix")
mt.genes <- which(rownames(data) %in% grep(pattern = "^mt-", x=rownames(data), value=TRUE))
percent.mito <- Matrix::colSums(data[mt.genes, ])/Matrix::colSums(data)
rp.genes <- which(rownames(data) %in% grep(pattern = "^Rp", x=rownames(data), value = TRUE))
percent.ribo <- Matrix::colSums(data[rp.genes, ])/Matrix::colSums(data)
data_nomt <- data[-mt.genes,]
data_norb <- data_nomt[-rp.genes,]

data <- CreateSeuratObject(counts=data_norb, min.cells = 3, min.features = 200, project = "10X_control_leukocytes", assay = "RNA")
data
data <- AddMetaData(object= data, metadata= percent.mito, col.name= "percent.mito")
data <- AddMetaData(object= data, metadata= percent.ribo, col.name= "percent.ribo")
length(data@meta.data$orig.ident)
data@meta.data$group <- rep("24hrs", times=3826)
data@meta.data$collection <- rep("2018", times=3826)
data@meta.data$dataset <- rep("CD45+ immune cells", times=3826)
data@meta.data$group2 <- rep("24hrs", times=3826)

VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol=4)

data <- NormalizeData(object=data, normalization.method = "LogNormalize", scale.factor=10000)
data <- FindVariableFeatures(object = data, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
head(x= HVFInfo(object = data))
data <- ScaleData(data)

data <- RunPCA(data, npcs = 20)
data <- RunUMAP(data, reduction="pca", dims=1:20)
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution=0.2)
head(Idents(data), 5)

DimPlot(data, pt.size = 0.5, label.size = 5, label=FALSE)
FeaturePlot(data, reduction = "umap", features = macrophages, ncol=5, pt.size = 0.5, order=TRUE, min.cutoff = 0)
#save RDS
saveRDS(data, file = "2018_T24_repairing_leukocytes.rds")

#=================================================================================================================================================================================================================
#read in seurat objects made with the above code for downstream integration and analysis
cycling_1 <- readRDS(file = "2018_cycling_leukocytes.rds")
cycling_2 <- readRDS(file = "2020_cycling_leukocytes.rds")
repairing_24hrs <- readRDS(file = "2018_T24_repairing_leukocytes.rds")
repairing_48hrs <- readRDS(file = "2020_T48_repairing_leukocytes.rds")

VlnPlot(object = repairing_48hrs, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), ncol=4, pt.size = 0)
DimPlot(repairing_48hrs, pt.size = 0.5, label.size = 5, label=TRUE)
FeaturePlot(cycling_2, reduction = "umap", features = "Krt18", ncol=4, pt.size = 0.5, order=TRUE, min.cutoff = 0)

#filtering seurat objects before integration
cycling_1 <- subset(x=cycling_1, subset= nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 0.1)
cycling_2 <- subset(x=cycling_2, subset= nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 0.1)
repairing_24hrs <- subset(x=repairing_24hrs, subset= nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 0.1)
repairing_48hrs <- subset(x=repairing_48hrs, subset= nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 0.1)

repairing_48hrs <- RenameIdents(repairing_48hrs, "0"="Neuts", "1"="Mono/Mac", "2"="Neuts", "3"="NK", "4"="cDC1/2", "5"="T", "6"="remove", 
                                "7"="Mono/Mac", "8"="B", "9"="remove", "10"="remove")
repairing_48hrs$celltype_1 <- Idents(repairing_48hrs)
levels(repairing_48hrs) <- c("Neuts", "Mono/Mac", "cDC1/2", "NK", "T", "B", "remove" )
repairing_48hrs <- subset(repairing_48hrs, idents = c("Neuts", "Mono/Mac", "cDC1/2", "NK", "T", "B"))

DimPlot(repairing_48hrs, pt.size = 0.5, label.size = 5, label=TRUE)
saveRDS(cycling_2, file = "2020_cycling_leukocytes_filtered.rds")

#integrate datasets together
anchors <- FindIntegrationAnchors(object.list=c(cycling_1, cycling_2, repairing_24hrs, repairing_48hrs), dims=1:20)
all_genes <- lapply(c(cycling_1, cycling_2, repairing_24rs, repairing_48hrs), row.names)  %>% Reduce(intersect, .)

leukocytes_combined <- IntegrateData(anchorset = anchors, new.assay.name = "integrated", features.to.integrate = all_genes)
DefaultAssay(leukocytes_combined) <- "integrated"
#leukocytes_combined <- subset(x=leukocytes_combined, subset= nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mito < 0.1 & doublet_score < 0.4)

leukocytes_combined <- ScaleData(leukocytes_combined)
leukocytes_combined <- RunPCA(leukocytes_combined, npcs = 20)
leukocytes_combined <- RunUMAP(leukocytes_combined, reduction="pca", dims=1:20)
leukocytes_combined <- FindNeighbors(leukocytes_combined, dims = 1:20)
leukocytes_combined <- FindClusters(leukocytes_combined, resolution=0.5)
head(Idents(leukocytes_combined), 5)

DimPlot(leukocytes_combined, pt.size = 0.5, label.size = 5, label=TRUE, group.by = "celltype_2")
FeaturePlot(leukocytes_combined, reduction = "umap", features = cDC2, ncol=4, pt.size = 0.5, order=TRUE, min.cutoff = 0)
leukocytes_combined$group <- factor(x=leukocytes_combined$group, levels = c("Cycling", "24hrs", "48hrs"))
FeatureScatter(leukocytes_combined, feature1 = "nFeature_RNA", feature2 = "Arg1")

#use gene expression to identify clusters
neutrophils <- c("S100a8","S100a9", "Cxcr2")
macrophages <- c("Csf1r", "Cx3cr1", "Fcgr1", "Adgre1", "Zeb2", "Ccr2", "Ly6c2", "H2-Ab1", "Cd68", "Itgam", "Cd14", "Lyz2", "Arg1", "Cd93", "Vegfa")
macrophages2 <- c("Csf1r", "Cx3cr1", "Fcgr1", "Adgre1", "Zeb2", "Ly6c2","Cd68", "Itgam")
NK <- c("Nkg7", "Gzma", "Gzmb")
cDC1 <- c("Xcr1", "Clec9a", "Btla", "Cadm1")
cDC2 <- c("Cd209a", "Clec10a")
Tc <- c("Cd3e", "Cd3d", "Cd3g")
Bc <- c("Cd79a", "Cd79b", "Ighm")
add <- c("Rpl10", "Rpl13")

FeaturePlot(leukocytes_combined, reduction = "umap", features = "Mrc1", ncol=1, pt.size = 0.5, order=TRUE, min.cutoff = 0)
VlnPlot(object = leukocytes_combined, features = macrophages2, ncol = 4, assay = "RNA", pt.size = 0)

DotPlot(leukocytes_combined, assay = "RNA", features = c(neutrophils, macrophages, cDC1, cDC2, NK, Tc, Bc), cols = c("lavender", "darkblue"), dot.min = 0, col.min = 0, dot.scale = 6, scale.by = "radius") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text((size=10))) + ggtitle("Canonical immune cell subtype markers")


#rename and reorder the clusters
leukocytes_combined <- RenameIdents(leukocytes_combined, "0"="Neuts_1", "1"="Neuts_3", "2"="Neuts_2", "3"="NK_1",
                                    "4"="Mac_1", "5"="Neuts_1", "6"="Mac_2", "7"="cDC2", "8"= "NK_2", "9"="Mac_3", 
                                    "10"="cDC3", "11"="Neuts_1", "12"="T", "13"="Mono", "14"="cDC1", "15"="B")

levels(leukocytes_combined) <- c("Neuts_1", "Neuts_2", "Neuts_3", "Mono", "Mac_1", "Mac_2", "Mac_3", 
                                 "cDC1", "cDC2", "cDC3", "NK_1", "NK_2", "T", "B")

leukocytes_combined$celltype_2 <- Idents(leukocytes_combined)
DimPlot(leukocytes_combined, reduction = "umap", pt.size = 0.5, label = FALSE, label.size = 6, split.by = "group")

saveRDS(leukocytes_combined, file = "2021_leukocytes_combined.Rds")
leukocytes_combined <- readRDS(file= "2021_leukocytes_combined.Rds")
#do doublet analysis
#=================================================================================================================================================================================================================
#DEgene analysis
#find markers for every cluster compared to all remaining cells, report only the positive ones
data.markers <- FindAllMarkers(leukocytes_combined, only.pos=TRUE, min.pct = 0.1, logfc.threshold = 0.5)
data.markers %>% group_by(cluster) %>% top_n(n=2, wt=avg_log2FC)
write.csv(data.markers, file= "leukocytes_cluster_markers.csv")

top20 <- data.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(leukocytes_combined, features = top20$gene)+ theme(axis.text.y = element_text(size=7))

#================================================================================================================================================================================================================
#subset to cell clusters of interest
#subset with cDC cells
leukocytes_subset <- subset(leukocytes_combined, idents = c("Mono", "Mac_1", "Mac_2", "Mac_3", "cDC1", "cDC2", "cDC3"))
leukocytes_subset <- subset(leukocytes_subset, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8"))

leukocytes_subset <- ScaleData(leukocytes_subset)
leukocytes_subset <- RunPCA(leukocytes_subset, npcs = 10)
leukocytes_subset <- RunUMAP(leukocytes_subset, reduction="pca", dims=1:10, min.dist = 0.2)
leukocytes_subset <- FindNeighbors(leukocytes_subset, dims = 1:10)
leukocytes_subset <- FindClusters(leukocytes_subset, resolution=0.3)
head(Idents(leukocytes_subset), 5)

DimPlot(leukocytes_subset, pt.size =1, label.size = 8, label=TRUE)
FeaturePlot(leukocytes_subset, reduction = "umap", features = "Casp3",split.by = "group", ncol=4, pt.size = 0.5, order=TRUE, min.cutoff = 0, cols = c("lavender", "darkblue"))
VlnPlot(object = leukocytes_subset, features = c("Csf1r","Zeb2", "S100a8", "S100a9"), ncol = 4, assay = "RNA", pt.size = 0)

interest <- c("Csf1r", "Zeb2", "Fcgr1", "Cd68", "Itgam", "Itgax", "Ccr2", "Cx3cr1", "Ly6c2", 
              "H2-Ab1","Lyz2","C1qa", "Cd72", "Adgre1", "Mertk",   "Cd14", "Cd93", "Arg1", "Vegfa", "Ccl2")

gating <- c("Ptprc", "Csf1r", "Itgam", "")
FeatureScatter(leukocytes_subset, feature1 = "Sell", feature2 = "Il1r2", group.by = "celltype_3", slot = "data")

#rename the clusters
leukocytes_subset <- RenameIdents(leukocytes_subset, "0"="Mac1", "1"="cDC2", "2"="Mac3", "3"="cDC3", "4"="Mac2", "5"="Mono", "6"="cDC1", "7"="Mac5", "8"="Mac4")
levels(leukocytes_subset) <- c("Mono", "Mac1", "Mac2", "Mac3", "Mac4", "Mac5", "cDC1", "cDC2", "cDC3")
leukocytes_subset$celltype_3 <- Idents(leukocytes_subset)

saveRDS(leukocytes_subset, file = "2021_leukocytes_subset.Rds")
leukocytes_subset <- readRDS(file = "2021_leukocytes_subset.Rds")

data.markers <- FindAllMarkers(leukocytes_subset, only.pos=TRUE, min.pct = 0.5, logfc.threshold = 0.5)
data.markers %>% group_by(cluster) %>% top_n(n=2, wt=avg_log2FC)
write.csv(data.markers, file= "leukocytes_subset_markers.csv")

top20 <- data.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
DoHeatmap(leukocytes_subset, features = top20$gene)+ theme(axis.text.y = element_text(size=7))

genes <- unique(top20$gene)
DotPlot(leukocytes_combined, assay = "RNA", features = genes, cols = c("lavender", "darkblue"), dot.min = 0, col.min = 0, dot.scale = 5, scale.by = "radius") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text((size=10))) + ggtitle("Cell subtype markers")

FeatureScatter(leukocytes_subset, feature1 = "Cd14", feature2 = "Siglec1", slot = "data", group.by = "celltype_3")

mono <- data.markers[data.markers$cluster=="Mono", 7]

interest <- top20$gene
interest <- unique(interest)

DotPlot(leukocytes_subset, assay = "RNA", features = interest, cols = c("lavender", "darkblue"), dot.min = 0, col.min = 0, dot.scale = 6, scale.by = "radius") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10), axis.text.y = element_text((size=10))) + ggtitle("Mono/Mac subtype markers")

#cell numbers and proportions
write.csv(table(Idents(leukocytes_subset), leukocytes_subset$group), file = "cell_numbers_subset.csv")
write.csv(prop.table(table(Idents(leukocytes_subset), leukocytes_subset$group), margin = 2), file = "cell_props_subset.csv")

#================================================================================================================================================================================================================
#analysis of surface receptors for each cluster
surface <- read.csv(file="GO gene lists of interest/mouse_cell_surface_receptors_CSPA.csv")

data.markers$surface <- surface$ENTREZ.gene.symbol[match(data.markers$gene, surface$ENTREZ.gene.symbol)]
data.markers.surface <- data.markers[complete.cases(data.markers), ]
data.markers.surface <- data.markers.surface[order(data.markers.surface$pct.2, decreasing= FALSE), ]

write.csv(data.markers.surface, file = "cell_surface_proteins.csv")
data.markers.surface <- read.csv(file = "leukocyte_analysis/subset_cDCs_mono_mac/cell_surface_proteins.csv")
surface_hits <- unique(data.markers.surface$gene)

DotPlot(leukocytes_subset, assay = "RNA", features =surface_hits, cols = c("lavender", "darkblue"), dot.min = 0, col.min = 0, dot.scale = 5, scale.by = "radius") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10), axis.text.y = element_text((size=10))) + ggtitle("Surface proteins")
DoHeatmap(leukocytes_subset, features = surface_hits)+ theme(axis.text.y = element_text(size=7))

mono <- data.markers.surface[data.markers.surface$cluster == "Mono", 7]
mac1 <- data.markers.surface[data.markers.surface$cluster == "Mac1", 7]
mac2 <- data.markers.surface[data.markers.surface$cluster == "Mac2", 7]
mac3 <- data.markers.surface[data.markers.surface$cluster == "Mac3", 7]
mac4 <- data.markers.surface[data.markers.surface$cluster == "Mac4", 7]
mac5 <- data.markers.surface[data.markers.surface$cluster == "Mac5", 7]
cDC1 <- data.markers.surface[data.markers.surface$cluster == "cDC1", 7]
cDC2 <- data.markers.surface[data.markers.surface$cluster == "cDC2", 7]
cDC3 <- data.markers.surface[data.markers.surface$cluster == "cDC3", 7]

chosen <- c("Pdpn", "Thbs1", "Cd14", "Thbd",
            "Fcgr1", "Siglec1", "C3ar1", "Msr1" ,
            "Cd72", "C1qa", "Tgfbr1", "Abca1", 
            "Sell", "Emb", "Il17ra", "Il6ra", 
            "Smpdl3b", "Slc12a2", "Megf9", "Ace", 
            "Il1r2", "Csf1", "Cd33", "Il1rap",
            "Hepacam2", "Itgae", "Anpep", "Cadm1", 
            "Olfm1", "Itgb7", "Itpr1", "Alcam", 
            "Il7r", "Cd200", "Ccr7", "Tspan3")

FeaturePlot(leukocytes_subset, reduction = "umap", features = "Jag1", split.by = "group", pt.size = 0.5, order=TRUE, min.cutoff = 1)
VlnPlot(object = leukocytes_subset, features = "Cd93", assay = "RNA", pt.size = 0)

#calculate cluster averages
cluster.averages <- AverageExpression(leukocytes_subset)
cluster.averages <- as.data.frame(cluster.averages$RNA)

surface1 <- cluster.averages[rownames(cluster.averages) %in% data.markers.surface$gene, ]

#generate clustered heatmap of average expression
subset_matrix <- surface1 %>% t() %>% scale() %>% t()

col <- colorRampPalette(c("white", "lavender", "darkblue"))(20)
heatmap(subset_matrix, col=col, cexRow = 0.5, cexCol = 1, main="Surface proteins") 
out <- heatmap(subset_matrix, col=col, cexRow = 0.5, cexCol = 1, main="Surface proteins") 

#reorder to match heatmap
h_order <- rev(rownames(subset_matrix)[out$rowInd])
reorder <- match(h_order, rownames(surface1))
surface1_reordered <- surface1[reorder, ]

#write csv file
write.csv(surface1_reordered, file = "all_celltypes_surface_proteins_ordered.csv")

#analysis of transcription factors for each cluster
tfs <- read.csv(file="GO gene lists of interest/mouse_TF_list.csv")
data.markers$tf <- tfs$Symbol[match(data.markers$gene, tfs$Symbol)]
data.markers.tf <- data.markers[complete.cases(data.markers), ]

tf_hits <- unique(data.markers.tf$gene)

DotPlot(leukocytes_subset, assay = "RNA", features =tf_hits, cols = c("lavender", "darkblue"), dot.min = 0, col.min = 0, dot.scale = 5, scale.by = "radius") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text((size=10))) + ggtitle("Transcription factors")

mono <- data.markers[data.markers$cluster=="Mono", 7]
#================================================================================================================================================================================================================

#Gene ontology for differentially expressed genes(clusterProfiler)
library(clusterProfiler)
library("org.Mm.eg.db")

#input from SCseq: degenes per cluster: data.markers <- read.csv(file= "all_cluster_markers.csv")
#add ensembl IDs to data frame
gene_ensembl <- read.table("RAW_FILES/2018_control_leukocytes/outs/filtered_feature_bc_matrix/features.tsv.gz")
data.markers$ensembl <- gene_ensembl$V1[match(data.markers$gene, gene_ensembl$V2)]
data.markers$celltype <- data.markers$cluster
head(data.markers)

#formatting gene lists to match input for clusterProfiler using data.markers file and GeneFormatClusterProfiler function
Mono_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "Mono")
Mac1_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "Mac1")
Mac2_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "Mac2")
Mac3_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "Mac3")
Mac4_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "Mac4")
Mac5_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "Mac5")
cdc1_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "cDC1")
cdc2_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "cDC2")
cdc3_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "cDC3")

Mono_genes <- names(Mono_genelist)
Mac1_genes <- names(Mac1_genelist)
Mac2_genes <- names(Mac2_genelist)
Mac3_genes <- names(Mac3_genelist)
Mac4_genes <- names(Mac4_genelist)
Mac5_genes <- names(Mac5_genelist)
cdc1_genes <- names(cdc1_genelist)
cdc2_genes <- names(cdc2_genelist)
cdc3_genes <- names(cdc3_genelist)

#individual cluster analysis using clusterProfiler- GOBP, GOCC, GOMF 
GOBP <- enrichGO(cdc3_genes, ont="BP", OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
dotplot(GOBP, showCategory=20) + ggtitle("GOBP:cDC3")
head(summary(GOBP))

#gene concept network
GOBP_network <- setReadable(GOBP, "org.Mm.eg.db", "ENSEMBL")
cnetplot(GOBP_network, showCategory = 10, categorySize= "pvalue", foldChange = cdc3_genelist, colorEdge=TRUE, node_label="gene")
#heatmap classification
heatplot(GOBP_network, showCategory= 20, foldChange = cdc3_genelist)+
  ggplot2::scale_fill_gradient(low = "blue", high = "red")

#list of lists; "V"=V_genes, "P"=P_genes, "F1"=F1_genes, "F2"=F2_genes, "F3"=F3_genes, "E"=E_genes
all_genelist2 <- list("Mono"= Mono_genes, "Mac1"=Mac1_genes, "Mac2"=Mac2_genes, "Mac3"=Mac3_genes, "Mac4"=Mac4_genes, "Mac5"=Mac5_genes)
#running compareCluster
All_GOBP2 <- compareCluster(geneClusters = all_genelist2, fun = "enrichGO", ont="BP", OrgDb= "org.Mm.eg.db", keyType = "ENSEMBL")
head(summary(All_GOBP))

All_GOBP2_2 <- simplify(All_GOBP2, cutoff=0.7, by="p.adjust", select_fun=min)

dotplot(All_GOBP2_2, showCategory=15)

#trajectory analysis on subset data
library()
library(htmlwidgets)
library(SingleCellExperiment)

#BUILDING THE CDS FROM A SEURAT OBJECT: 
data <- leukocytes_subset
#gene annotations
gene_annotation <- as.data.frame(rownames(data@reductions[["pca"]]@feature.loadings), row.names = rownames(data@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
head(gene_annotation)
#cell information including group, seurat clusters and celltypes
cell_metadata <- as.data.frame(data@assays[["RNA"]]@counts@Dimnames[[2]], row.names=data@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
head(cell_metadata)
clusters <- as.data.frame(data$seurat_clusters)
colnames(clusters) <- "seurat_clusters"
cell_metadata <- cbind(cell_metadata, clusters)
names <- as.data.frame(data$celltype_3)
colnames(names) <- "celltype"
cell_metadata <- cbind(cell_metadata, names)
group <- as.data.frame(data$group)
colnames(group) <- "group"
cell_metadata <- cbind(cell_metadata, group)
cellgroup <- as.data.frame(data$cellgroup)
colnames(cellgroup) <- "cell_group"
cell_metadata <- cbind(cell_metadata, cellgroup)
head(cell_metadata)

#counts sparse matrix
new_matrix <- data@assays[["RNA"]]@counts
new_matrix <- new_matrix[rownames(data@reductions[["pca"]]@feature.loadings),]
expression_matrix <- new_matrix

#constructing the basic cds object
cds <- new_cell_data_set(expression_matrix, 
                                    cell_metadata = cell_metadata, 
                                    gene_metadata = gene_annotation)

#PRE-PROCESSING AND BUILDING TRAJECTORY
#Assign UMAP co-ordinates from seurat object
umap <- data@reductions[["umap"]]@cell.embeddings
colnames(umap) <- c("V1", "V2")
cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- umap

#cluster the cells
cds <- preprocess_cds(cds, num_dim = 20)
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "partition", cell_size = 1)

#learn the trajectory graph
cds  <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE, learn_graph_control = list(ncenter=300, rann.k=50))
plot_cells(cds, label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 4.5, cell_size = 0.75, color_cells_by = "celltype")

#ordering the cells along the pseudotime and visualise
cds  <- order_cells(cds , reduction_method = "UMAP")
plot_cells(cds , color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = TRUE, label_branch_points = TRUE, graph_label_size = 4, cell_size = 0.5)

saveRDS(cds, file = "2021__leukoctes_subset.Rds")

#================================================================================================================================================================
#specific degene comparisons: each celltype between timepoints
mono_cons <- FindConservedMarkers(leukocytes_subset, ident.1 = "Mono", grouping.var = "group")

#monocytes
mono_cells <- subset(leukocytes_subset, idents= "Mono")

Idents(mono_cells) <- "group"
avg.mono.cells <- as.data.frame(log1p(AverageExpression(mono_cells, verbose = FALSE)$RNA))
avg.mono.cells$gene <- rownames(avg.mono.cells)

genes <- avg.mono.cells$gene
A <- avg.mono.cells[avg.mono.cells$`24hrs` > ((avg.mono.cells$`48hrs`) + 1.5), 4]
B <- avg.mono.cells[avg.mono.cells$`48hrs` > ((avg.mono.cells$`24hrs`) + 1.5), 4]
C <- c(A, B)

ggplot(avg.mono.cells, aes(`48hrs`, `24hrs`)) + geom_point(col=ifelse(genes %in% C, "red", "grey"), size=ifelse(genes %in% C, 1, 0.5)) +ggtitle("Mono: average gene expression")+
  geom_text_repel(colour = "darkblue", size=3, aes(label=ifelse(genes %in% C, as.character(genes), '')))

#macrophages
mac1_cells <- subset(leukocytes_subset2, idents= "Mac_1")
Idents(mac1_cells) <- "group"
mac2_cells <- subset(leukocytes_subset2, idents= "Mac_2")
Idents(mac2_cells) <- "group"
mac3_cells <- subset(leukocytes_subset2, idents= "Mac_3")
Idents(mac3_cells) <- "group"

avg.mac1.cells <- as.data.frame(log1p(AverageExpression(mac1_cells, verbose = FALSE)$RNA))
avg.mac1.cells$gene <- rownames(avg.mac1.cells)

avg.mac2.cells <- as.data.frame(log1p(AverageExpression(mac2_cells, verbose = FALSE)$RNA))
avg.mac2.cells$gene <- rownames(avg.mac2.cells)

avg.mac3.cells <- as.data.frame(log1p(AverageExpression(mac3_cells, verbose = FALSE)$RNA))
avg.mac3.cells$gene <- rownames(avg.mac3.cells)

genes <- avg.mac3.cells$gene
A <- avg.mac3.cells[avg.mac3.cells$`24hrs` > ((avg.mac3.cells$`48hrs`) + 1.5), 4]
B <- avg.mac3.cells[avg.mac3.cells$`48hrs` > ((avg.mac3.cells$`24hrs`) + 1.5), 4]
C <- c(A, B)

ggplot(avg.mac3.cells, aes(`24hrs`, `48hrs`)) + geom_point(col=ifelse(genes %in% C, "red", "grey"), size=ifelse(genes %in% C, 1, 0.5)) +ggtitle("Mac3: average gene expression")+
  geom_text_repel(colour = "darkblue", size=3, aes(label=ifelse(genes %in% C, as.character(genes), '')))

#backplotting marker sets
mono1 <- c("Arg1", "Pdpn", "Mmp12", "Inhba")
mono2 <- c("Il2rb", "Nkg7", "Klra4", "Gzmb")
mac1 <- c("Fcgr1", "Ifit2", "Aif1", "Irf7")
mac2 <- c("Itga4", "Klf4", "Adgre4", "Gpr141")
mac3 <- c("Cd72", "Ms4a7", "C1qa", "H2-Ab1")
mac4 <- c("S100a9", "S100a8", "Il1b", "Cxcr2")

FeaturePlot(leukocytes_combined, reduction = "umap", features = "Mertk", ncol=4, pt.size = 0.5, order=TRUE, min.cutoff = 0)
FeatureScatter(leukocytes_subset, feature1 = "H2-Ab1", feature2 = "Ly6c2", slot = "data")

#pie charts
Idents(leukocytes_subset) <- leukocytes_subset$group
cycling <- subset(leukocytes_subset, idents="Cycling")
repair <- subset(leukocytes_subset, idents="24hrs")
resolve <- subset(leukocytes_subset, idents="48hrs")

Idents(resolve) <- resolve$celltype_3

A <- as.data.frame(prop.table(table(cycling$celltype_3)))
B <- as.data.frame(prop.table(table(repair$celltype_3)))
C <- as.data.frame(prop.table(table(resolve$celltype_3)))

colnames(C) <- c("group", "value")

bp <- ggplot(A, aes(x="Cycling", y=value, fill=group))+
  geom_bar(width=1, stat = "identity")
bp

pie <- bp + coord_polar("y", start = 0)
pie 

pie + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(label= percent(value)), size=4, position = position_stack(vjust=0.5))

pie + blank_theme +
  theme(axis.text.x=element_blank())

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    panel.border = element_blank(), 
    panel.grid = element_blank(), 
    axis.ticks = element_blank(), 
    plot.title = element_text(size=20, face="bold")
  )

#=============================================================================================
#FIGURE THIS OUT!!!! TRING TO CLUSTER ISAACS GENE LIST BY CELLTYPE AS WITH A HEATMAP AND THEN PRODUCE DOTPLOTS!
col <- colorRampPalette(c("white", "lavender", "darkblue"))(20)

surface <- read.csv(file = "GO gene lists of interest/mouse_cell_surface_receptors_CSPA.csv")
surface <- wounds[!duplicated(wounds$ENTREZ.gene.symbol), ]

data.markers.surface

data@meta.data

#calculate cluster averages
cluster.averages <- AverageExpression(leukocytes_subset)
cluster.averages <- as.data.frame(cluster.averages$RNA)

subset <- cluster.averages[rownames(cluster.averages) %in% data.markers.surface$gene, ]

#generate clustered heatmap of average expression
subset_matrix <- subset %>% t() %>% scale() %>% t()
subset_matrix[is.na(subset_matrix)] <- 0

heatmap(subset_matrix, col=col, cexRow = 0.5, cexCol = 2, main="Counts") 
out <- heatmap(subset_matrix, col=col, cexRow = 0.5, cexCol = 0.5, main="WT versus 5ar; decidualised") 

#reorder to match heatmap
h_order <- rev(rownames(subset_matrix)[out$rowInd])
reorder <- match(h_order, rownames(subset))
subset_reordered <- subset[reorder, ]

#write csv file
write.csv(subset_reordered, file = "leukocyte-surface_proteins_ordered.csv")

met <- c("Bmp4", "Cemip", "Cited1", "Ctnnb1", "Ehf", "Elf3", "Fzd7", "Gata3", "Grem1", "Hmga2", "Irf6", "Klf4", "Pax2", "Pax8", "Smo", "Wnt4", "Wt1")
emt <- c("Notch1", "Snai2", "Sox9", "Twist1")


            