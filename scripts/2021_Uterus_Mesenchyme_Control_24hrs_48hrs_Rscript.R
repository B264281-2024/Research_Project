library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
library(UpSetR)
source("SCRIPTS/Utils.R")
theme_set(theme_cowplot())

#read in seurat objects
mes <- readRDS("Seurat Objects/2020_cycling_mes.Rds")
epi <- readRDS("Seurat Objects/2020_cycling_epi.Rds")
mes_24 <- readRDS("Seurat Objects/2020_T24_mes.Rds")
mes_48 <- readRDS("Seurat Objects/2020_T48_mes.Rds")
epi_48 <- readRDS("Seurat Objects/2020_T48_epi.Rds")

#visualise outliers
DimPlot(mes)
  VlnPlot(object = epi_48, features = c("nFeature_RNA", "percent.mito", "doublet_score", "Pdgfrb", "Epcam"), ncol = 5, pt.size = 0.1)
FeatureScatter(object = mes_24, feature1 = "nFeature_RNA", feature2 = "percent.mito")
FeaturePlot(epi_T48, pt.size = 1, features = c("Krt5", "Krt14"), cols = c("lightgoldenrodyellow", "darkred"), min.cutoff = 0, label = TRUE, order = TRUE)

#subset if only interested in specific clusters
#mes
mes <- subset(mes, idents = c("V", "P", "F1", "F2", "F3"))
mes <- subset(mes, subset = (doublet_score < 0.5 & percent.mito < 0.1))
#mes_24
mes_24 <- subset(mes_24, idents = c("V", "P", "F1", "F2", "F3", "New1"))
mes_24 <- subset(mes_24, subset = (doublet_score < 0.5 & percent.mito < 0.2))
#mes_48
mes_48 <- subset(mes_48, subset = (doublet_score < 0.5 & percent.mito < 0.2))
#epi
epi <- subset(epi, subset = (doublet_score < 0.5 & Epcam > 2))
#epi_48
epi_48 <- subset(epi_T48, subset = (doublet_score < 0.5 & percent.mito < 0.1 & Epcam > 2))

ncol(mes_48)
mes_48@meta.data$timepoint <- rep("48hrs", times=6507)

#save filtered objects
saveRDS(mes_48, file = "Seurat Objects/2020_T48_mes_filtered.Rds")

##ANALYSIS OF MESENCHYMAL DATASETS ALONE
#integrate datasets together
anchors <- FindIntegrationAnchors(object.list=c(mes, mes_24, mes_48), dims=1:20)
all_genes <- lapply(c(mes, mes_24, mes_48), row.names)  %>% Reduce(intersect, .)

mes_combined <- IntegrateData(anchorset = anchors, new.assay.name = "integrated", features.to.integrate = all_genes)
mes_combined <- ScaleData(mes_combined)
mes_combined <- RunPCA(mes_combined, npcs = 30)
mes_combined <- RunUMAP(mes_combined, reduction="pca", dims=1:20, min.dist = 0.4)

mes_combined <- FindNeighbors(mes_combined, dims = 1:20)
mes_combined <- FindClusters(mes_combined, resolution=0.25)
head(Idents(mes_combined), 5)

mes_combined <- subset(mes_combined, idents = c("0", "1", "2", "3", "4", "5", "6"))

DimPlot(mes_combined, pt.size = 0.5, label.size = 5, label=TRUE, split.by = "timepoint")
mes_combined$timepoint <- factor(x=mes_combined$timepoint, levels = c("cycling", "24hrs", "48hrs"))

FeaturePlot(mes_combined, reduction ="umap", split.by = "timepoint", features= c("Pdgfrb", "Pdgfra", "Mcam", "Cspg4", "Hif1a"), min.cutoff = 0.25, ncol=5, pt.size=1, order=TRUE, cols=c("lavender", "blue"))
VlnPlot(mes_combined, features= c("Pdgfrb", "Pdgfra", "Mcam", "Cspg4", "Hif1a"),  assay="RNA", ncol=5, pt.size=0)
FeatureScatter(mes_combined, feature1 = "Pdgfrb", feature2 = "Krt18", group.by = "timepoint")
#=================================================================================================================================================================================================================
#running ChangeOrderAndIds to reorder and rename the clusters
mes_combined <- ChangeOrderAndIds(object = mes_combined,
                                     new_order= c("4", "3", "7","2", "9", "0", "1", "8", "5", "6"), 
                                     new_names = c("V", "P1", "P2", "F1", "F1", "F2", "F3", "F3", "F4", "F5"))
mes_combined@meta.data$celltype2 <- Idents(mes_combined)

#visualise
hue_pal()DimPlot(mes_combined, pt.size = 0.5, label.size = 6, label = TRUE)

#save seurat object
saveRDS(mes_combined, file = "Seurat Objects/2021_cycling_T24_T48_mesenchyme.Rds")
mes_combined <- readRDS(file = "Seurat Objects/2021_cycling_T24_T48_mesenchyme.Rds")

identificationV <- c("Acta2","Myh11", "Tagln", "Cnn1", "Kitl", "Pln", "Lmod1", "Sorbs2")
identificationP <- c("Rgs5", "Kcnj8", "Myo1b", "Abcc9", "Ednrb", "Vtn", "Nrp1", "Postn")
identificationF1 <- c("Spon2", "Aspg", "Dpep1", "Ngfr", "Angptl7", "Ifit1", "Ifit3", "H2-Q7")
identificationF2 <- c("Cxcl14", "Cdh11", "Wt1", "Rgs2", "Smoc2", "Wnt4", "Aldh1a2", "Bmp7")
identificationF3 <- c("Clec3b", "Col14a1", "Fap", "Cd55", "Cxcl16", "Mmp3", "Efemp1", "Vit")
stromal <- c("Pdgfra", "Cd34", "Fbln1", "Fbln2", "Col1a1", "Mfap4", "Mfap5")
pericyte <- c("Mcam", "Cspg4", "Rgs4", "Angpt2", "Ngf")
epithelial <- c("Cdh1", "Epcam", "Muc1", "Muc4", "Krt7", "Alcam", "Fxyd3", "Fxyd4", "Cd9")
canonical <- c("Vim", "Des", "Thy1", "Pdgfrb", "Col1a1", "Mcam", "Cspg4", "Rgs4","Rgs5", "Kcnj8", "Pdgfra", "Fbln2", "Mfap4", "Mfap5", "Cd34", 
               "Krt8", "Krt18", "Krt19", "Epcam", "Clu", "Sdc4")

epithelial2 <- c("Alcam","Car2", "Cd9", "Cd24a", "Cdh1", "Cldn4", "Clu", "Dsp","Epcam", "Folr1", "Fxyd3", "Fxyd4", 
                  "Krt4", "Krt5", "Krt7","Krt8", "Krt9", "Krt13", "Krt14", "Krt15", "Krt18",  "Krt19", 
                 "Lcn2", "Ly6a", "Muc1", "Muc4","Ocln","Olfm1","Plet1","Prap1", "Prom1", "Sdc4", "Sprr2a", "Sprr2f","Trpv6", "Wfdc2" )
epithelial3 <- c("Alcam","Car2", "Cd9", "Cd24a", "Cdh1", "Cldn4", "Clu", "Dsp","Epcam", "Fxyd3", "Fxyd4", 
                  "Krt5", "Krt7","Krt8", "Krt9", "Krt14", "Krt15", "Krt18",  "Krt19", 
                 "Lcn2", "Muc1", "Muc4","Olfm1","Plet1","Prap1", "Prom1", "Sdc4", "Sprr2a", "Sprr2f","Trpv6", "Wfdc2" )
mesenchymal <- c("Cd34","Cdh11","Col1a1", "Col3a1", "Col6a4", "Des","Dio2","Fbln1", "Fbln2", "Mfap4", "Mfap5", "Pdgfra", "Pdgfrb", "Thy1", "Vim")
MET <- c("Bmp4", "Cemip", "Cited1", "Ctnnb1", "Dspp", "Ehf", "Elf3", "Foxa1", "Fzd7", "Gata3","Gdnf","Grem1", "Grhl1",
         "Grhl2", "Grhl3", "Irf6", "Klf4","Lif", "Ovol2","Pax2", "Pax8", "Sall1", "Smo", "Tcf15","Trp63", 
         "Wnt4", "Wnt9b", "Wt1", "Ybx2", "Znf165")
stem <- c("Cd34", "Krt14", "Krt15", "Nestin", "Trp63", "Itgb1", "Sox9", "Lgr5", "Itga6", "S100a6")
decidualisation <- c("Ash1l", "Cdh1", "Cited2", "Ctsb", "Ctsl", "Dedd", "Epor", "Ghrl", "Ghsr", "Il11ra1", "Junb", "Kiss1" ,"Lif", "Mapk1", "Mapk3", "Pla2g4a", "Ptgis", "Ptgs2", "Ptn", "Tcf23", "Tppp3")
M_E <- c("Vim", "Pdgfrb", "Thy1", "Cd34",  "Pdgfra", "Cdh11", "Fbln1",  "Mfap4", "Krt8", "Krt18", "Cd24a", "Alcam", "Epcam", "Sdc4", "Plet1", "Clu")
F4_specific <- c("Mmp12", "Krt18", "Mt1", "Mt2", "Il23a", "Cyp26b1", "Tnfrsf11b", "Cemip", "grem2", "Grem1", "Hmga2", 
                 "Edn1", "Tnfrsf23", "Fam84a", "Nrg1", "Tnf", "Steap4", "Coch", "Clu", "Pdpn", "Ald1a3", "Sprr1a", 
                 "Smox", "Gsr")

wounds <- read.csv(file = "GO gene lists of interest/EMT_GO_term_summary_20210621_061710.csv")
wounds <- wounds[!duplicated(wounds$Symbol), ]
table(wounds$Annotated.Term)

chosen <- intersect(wounds$Symbol, data.markers$gene)

DotPlot(mes_combined, assay = "RNA", features = wounds$Symbol, cols = c("lavender", "darkblue"), dot.min = 0, col.min = 0, dot.scale = 6, scale.by = "radius") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text((size=10))) + ggtitle("Epithelial to mesenchymal transition")

FeaturePlot(mes_combined, reduction ="umap", features= "Thy1", min.cutoff = 1, ncol=1, pt.size=0.25, order=TRUE, cols=c("lavender", "darkblue"))
VlnPlot(mes_combined, assay = "RNA", features= c("Thy1"), ncol=1, pt.size=0)
StackedVlnPlot(mes_combined, assay = "RNA", features= c("Pdgfrb", "Vim", "Pdgfra", "Cdh11", "Krt18", "Sdc4", "Cd24a", "Epcam"))
#degene analysis
data.markers <- DEgeneAnalysis(object=mes_combined)
data.markers <- read.csv(file= "exploratory_analysis_timepoint_mes_epi/all_cluster_markers.csv")

#if reading back in use: data.markers <- read.csv(file= "all_cluster_markers.csv")

###PLOTS: visualisation of DEgene results using gene sets of interest (see gene lists in .csv file)
#heatmap of DEgenes per cluster
top20 <- data.markers %>% group_by(cluster) %>% top_n(n=15, wt=avg_log2FC)
F4 <- data.markers[data.markers$cluster=="F4", ]
F4 <- F4[order(F4$avg_log2FC, decreasing = TRUE), 7]

DoHeatmap(mes_combined, features = F4) + 
  NoLegend()+
  theme(axis.text.y = element_text(size=6))

chosen <- c("Mustn1", "Myh11", "Lmod1",
            "Kcnj8", "Vtn", "Abcc9",
            "Serpine1", "Arhgdib", "Klk8",
            "Dpep1", "Ngfr", "Sema3c", 
            "Cxcl14", "Smoc2", "Wt1", 
            "Lum", "Efemp1", "Col3a1", 
            "Krt18", "Mt2", "Pmaip1", 
            "Col7a1", "Scube1", "Neat1")

DotPlot(mes_combined, assay = "RNA", col.min = 0.5, features = F4[256:340], dot.scale = 4, cols = c("lavender", "darkblue"), dot.min = 0.1, scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text((size=120)))

VlnPlot(mes_combined, features = epithelial, ncol = 6, pt.size = 0, assay="RNA")
FeaturePlot(mes_combined, reduction ="umap", features= chosen, min.cutoff = 1, ncol=6, pt.size=0.5, order=TRUE, cols=c("lavender", "darkblue"))
FeaturePlot(mes_combined, reduction ="umap", features= c("Egfr", "Pdgfra", "Bcl2"), min.cutoff = 0, ncol=1, pt.size=1, order=TRUE, cols=c("lavender", "darkblue"))


#Gene ontology for differentially expressed genes(clusterProfiler)
library(clusterProfiler)
library("org.Mm.eg.db")

#input from SCseq: degenes per cluster: data.markers <- read.csv(file= "all_cluster_markers.csv")
#add ensembl IDs to data frame
gene_ensembl <- read.table("RAW_FILES/T24_control_mesenchyme/outs/filtered_feature_bc_matrix/features.tsv.gz")
data.markers$ensembl <- gene_ensembl$V1[match(data.markers$gene, gene_ensembl$V2)]
data.markers$celltype <- data.markers$cluster
head(data.markers)

#adding ensembl IDs to a gene list
#VP_genes$gene <- rownames(VP_genes)
#VP_genes$ensembl <- gene_ensembl$V1[match(VP_genes$gene, gene_ensembl$V2)]

#formatting gene lists to match input for clusterProfiler using data.markers file and GeneFormatClusterProfiler function
V_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "V")
P1_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "P1")
P2_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "P2")
F1_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "F1")
F2_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "F2")
F3_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "F3")
F4_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "F4")
F5_genelist <- GeneFormatClusterProfiler(object = data.markers, cluster= "F5")

V_genes <- names(V_genelist)
P1_genes <- names(P1_genelist)
P2_genes <- names(P2_genelist)
F1_genes <- names(F1_genelist)
F2_genes <- names(F2_genelist)
F3_genes <- names(F3_genelist)
F4_genes <- names(F4_genelist)
F5_genes <- names(F5_genelist)

#individual cluster analysis using clusterProfiler- GOBP, GOCC, GOMF 
GOBP_F4 <- enrichGO(F4_genes, ont="BP", OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GOBP_F4_2 <- gofilter(GOBP_F4, level = 4)
dotplot(GOBP_F4_2, showCategory=20)
head(summary(GOBP_F4))

#gene concept network
GOBP_F4_network <- setReadable(GOBP_F4_2, "org.Mm.eg.db", "ENSEMBL")
cnetplot(GOBP_F4_network, showCategory = 20, categorySize= "pvalue", foldChange = F4_genelist)
cnetplot(GOBP_F4_network, showCategory = 5, categorySize= "pvalue", foldChange = F4_genelist, colorEdge=TRUE, node_label="gene")
#heatmap classification
heatplot(GOBP_F4_network, showCategory= 20, foldChange = F4_genelist)+
  ggplot2::scale_fill_gradient(low = "blue", high = "red")

#list of lists; "V"=V_genes, "P"=P_genes, "F1"=F1_genes, "F2"=F2_genes, "F3"=F3_genes, "E"=E_genes
all_genelist <- list("V"=V_genes, "P1"=P1_genes, "P2"=P2_genes, "F1"=F1_genes, "F2"=F2_genes, "F3"=F3_genes, "F4"=F4_genes, "F5"=F5_genes)
all_genelist2 <- list("F1"=F1_genes, "F2"=F2_genes, "F3"=F3_genes, "F4"=F4_genes, "F5"=F5_genes)

#running compareCluster
All_GOBP2 <- compareCluster(geneClusters = all_genelist2, fun = "enrichGO", ont="BP", OrgDb= "org.Mm.eg.db", keyType = "ENSEMBL")
head(summary(All_GOBP))

All_GOBP2_2 <- simplify(All_GOBP2, cutoff=0.7, by="p.adjust", select_fun=min)
All_GOBP2_3 <- gofilter(All_GOBP2, level = 4)

dotplot(All_GOBP2_3, showCategory=10)

##REACTOMEPA analysis
#Gene ontology for differentially expressed genes(ReactomePA)
library("ReactomePA")

#mapping entrez gene ids onto gene lists
gnames_entrez_F4 <- mapIds(org.Mm.eg.db, F4_genes,"ENTREZID", "ENSEMBL")

#KEGG enrichment for individual clusters
KEGG_F4 <- enrichKEGG(gnames_entrez_F4, organism = "mmu", pvalueCutoff = 0.05 )
dotplot(KEGG_F4, showCategory=20,  title= "KEGG: F4 cells")

#reactomePA for individual clusters
reactome_F4 <- enrichPathway(gene = gnames_entrez_F4, organism = "mouse", pvalueCutoff = 0.05, readable = TRUE)
head(as.data.frame(reactome_F4))

dotplot(reactome_F4, showCategory=20, title= "ReactomePA: New cells")
cnetplot(reactome_F4, categorySize="pvalue", foldChange = V_versus_P)

# ===========================================================================================================================================================================================================
#using upsetR for analysis of gene lists
library(UpSetR)
#gene list attributed to each cell cluster from standard data.markers file
fib1 <- data.markers[data.markers$cluster=="F1", 7]
fib2 <- data.markers[data.markers$cluster=="F2", 7]
fib3 <- data.markers[data.markers$cluster=="F3", 7]
fib4 <- data.markers[data.markers$cluster=="F4", 7]
fib5 <- data.markers[data.markers$cluster=="F5", 7]
per1 <- data.markers[data.markers$cluster=="P1", 7]
vas <- data.markers[data.markers$cluster=="V", 7]
per2 <- data.markers[data.markers$cluster=="P2", 7]

ListInput <- list(F1= fib1, F2=fib2, F3=fib3, F4=fib4, F5=fib5, P1=per1, P2=per2, V=vas)
upset(fromList(ListInput) , nsets = 8, order.by = "freq", 
      mainbar.y.label = "Number of DE genes per cell type",
      nintersects = 32,
      cutoff = 4,
      group.by = "sets",
      sets.x.label = "DE genes", 
      sets= c("V", "P1", "P2", "F1", "F2", "F3", "F4", "F5"), 
      keep.order = TRUE,
      sets.bar.color = c("coral1","gold2","olivedrab2", "seagreen3", "aquamarine2", "royalblue1", "mediumpurple1", "hotpink"), 
      point.size = 3,  line.size = 1)
#queries = list(list(query= intersects, params= list("New", "F2"), color="violetred", active=TRUE), list(query=intersects, params=list("New", "E"), color="blue", active=TRUE)),

#BASIC INTERACTION WITH SEURAT OBJECT AND GENERATING USEFUL MATRICES FOR DOWNSTREAM QUERIES
#tabulate cells by cluster ID, dataset or both
table(Idents(mes_combined))
table(mes_combined$group)
prop.table(table(Idents(mes_combined)))
table(Idents(mes_combined), mes_combined$timepoint)
prop.table(table(Idents(mes_combined), mes_combined$timepoint))

#gathering cells per cluster
V_cells <- WhichCells(mes_combined, idents = "V")

#extract the expression matirx for each cell cluster
F5_expression <- as.matrix(GetAssayData(mes_combined, slot = "data")[, WhichCells(mes_combined, idents = "F5")])

#calculating the average expression of all cells within a cluster
cluster.averages <- AverageExpression(mes_combined)
head(cluster.averages[["RNA"]], 100)
cluster_expression <- as.data.frame(cluster.averages$RNA)

V_average <- rownames(cluster_expression[cluster_expression$V > 1, ])
P1_average <- rownames(cluster_expression[cluster_expression$P1 > 1, ])
P2_average <- rownames(cluster_expression[cluster_expression$P2 > 1, ])
F1_average <- rownames(cluster_expression[cluster_expression$F1 > 1, ])
F2_average <- rownames(cluster_expression[cluster_expression$F2 > 1, ])
F3_average <- rownames(cluster_expression[cluster_expression$F3 > 1, ])
F4_average <- rownames(cluster_expression[cluster_expression$F4 > 1, ])
F5_average <- rownames(cluster_expression[cluster_expression$F5 > 1, ])

length(F1_average)
head(F1_average, 20)

ListInput <- list(F1= F1_average, F2=F2_average, F3=F3_average, F4=F4_average, F5=F5_average, P1=P1_average, P2=P2_average, V=V_average)
upset(fromList(ListInput) , nsets = 8, order.by = "freq", 
      mainbar.y.label = "Number of expressed genes",
      nintersects = 50,
      cutoff = 4,
      sets.x.label = "Gene expression > 1", 
      sets= c("V", "P1", "P2", "F1", "F2", "F3", "F4", "F5"), 
      keep.order = TRUE,
      sets.bar.color = c("coral1","gold2","olivedrab2", "seagreen3", "aquamarine2", "royalblue1", "mediumpurple1", "hotpink"), 
      point.size = 3,  line.size = 1 
      )

#finding unique genes to each cluster
to.remove <- F4_average %in% unique(c(V_average, P1_average, F1_average, F2_average, F3_average))
F4_unique <- as.list(F4_average[!to.remove])

DotPlot(mes_combined, assay = "RNA", col.min = 0, features = F4_unique, dot.scale = 4, cols = c("lavender", "darkblue"), dot.min = 0.1, scale.by = "size") + RotatedAxis()+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8), axis.text.y = element_text((size=120)))

#quick GO analysis of new gene lists
#add ensembl IDs to data frame
gene_ensembl <- read.table("RAW_FILES/T24_control_mesenchyme/outs/filtered_feature_bc_matrix/features.tsv.gz")
F4_unique_genes <- gene_ensembl$V1[match(F4_unique, gene_ensembl$V2)]

#individual cluster analysis using clusterProfiler- GOBP, GOCC, GOMF 
GOBP_F4_unique <- enrichGO(F4_unique_genes, ont="BP", OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
GOBP_F4_unique2 <- gofilter(GOBP_F4, level = 4)
dotplot(GOBP_F4_unique, showCategory=20)
head(summary(GOBP_F4_unique))

F4_genes

# barplots of expression and statistics
#subset based on krt18
GOI <- FetchData(mes_combined, vars = c("Krt18"))
krt18 <- cluster_expression["Krt18",]
