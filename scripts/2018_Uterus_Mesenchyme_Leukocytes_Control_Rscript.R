#Analysis using SingleCellSignalR to infer intercellular networks from singl-cell transcriptomics
#date:11/06/2020

#data format: table of read or UMI counts, one per column per individual cell and one row per gene
BiocManager::install ("Seurat")
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
library(UpSetR)
library(VennDiagram)
theme_set(theme_cowplot())
source("Utils.R")

#control mesenchyme and leukocytes
control.mesenchyme <- readRDS(file = "Seurat Objects/Control_Mesenchyme_seurat_object")
control.leukocytes <- readRDS(file = "Seurat Objects/Control_leukocytes_original")

control.mesenchyme@meta.data$type <- rep("mesenchyme", times=6275)
control.leukocytes@meta.data$type <- rep("leukocytes", times=1358)

#merge the two?
data.merged <- merge(control.mesenchyme, y=control.leukocytes, add.cell.ids = c("mesenchyme", "leukocytes"), project = "control M/L merge")
head(colnames(data.merged))
table(data.merged$orig.ident)
data.merged

#tidy up(remove Mt and Rp genes)
mt.genes <- which(rownames(data.merged) %in% grep(pattern = "^mt-", x=rownames(data.merged), value=TRUE))
rp.genes <- which(rownames(data.merged) %in% grep(pattern = "^Rp", x=rownames(data.merged), value = TRUE))
combined.data_nomt <- data.merged[-mt.genes,]
combined.data_norb <- combined.data_nomt[-rp.genes,]
data.merged <- combined.data_norb

data.merged <- ScaleData(data.merged)
data.merged <- FindVariableFeatures(object=data.merged)
data.merged <- RunPCA(data.merged, npcs=30, features = VariableFeatures(object = data.merged))
data.merged <- RunUMAP(data.merged, reduction= "pca", dims=1:20)
data.merged <- FindNeighbors(data.merged, dims = 1:20)
data.merged <- FindClusters(data.merged, resolution=0.1)
head(Idents(data.merged), 5)

data$celltype <- factor(data$celltype, levels= c("V", "P", "F1", "F2", "F3", "E", "Neuts", "Macs", "cDC1", "cDC2", "NK", "Tcells", "Bcells"))
DimPlot(data, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 6, group.by = "celltype")

saveRDS(data, file = "Control_mesenchyme_leukocytes_merged_Seurat_Object.rds")

#Using NicheNetR to perform ligand-receptor analysis
library(nichenetr)
library(tidyverse)

#check data
seuratObj <- readRDS(file = "Control_Repairing_Mesenchyme_Leukocytes_Merged_Seurat_Object.rds")
DimPlot(seuratObj, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 6, group.by = "celltype")
seuratObj@meta.data %>% head()

#visualise which cell types are present in the data
seuratObj@meta.data$celltype %>% table()
seuratObj@meta.data$group %>% table()

#read in NicheNet's ligand-target prior model, ligand-receptor network and weighted integrated networks
ligand_target_matrix <- readRDS(url("http://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
ligand_target_matrix[1:5, 1:5]

lr_network <- readRDS(url("http://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks <- readRDS(url("http://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from, to), by = c("from", "to"))
head(weighted_networks$lr_sig)
head(weighted_networks$gr)

#convert Nichenet network gene symbols form human to mouse based on one-to-one orthology
lr_network <- lr_network %>% mutate(from= convert_human_to_mouse_symbols(from), to= convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) <- ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) <- ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix <- ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr <- weighted_networks_lr %>% mutate(from= convert_human_to_mouse_symbols(from), to= convert_human_to_mouse_symbols(to)) %>% drop_na()
head(weighted_networks_lr)

##PERFORM THE NICHENET ANALYSIS
#apply NicheNet to predict which ligands expressed by all cells in the dataset are most likely to have induced the differential expression in CELLS OF INTEREST after condition change (ie during repair)
#NicheNet analysis pipeline
#1. Define a "sender/niche" cell population and a "receiver/target" cell population present in your expression data and determine which genes are expressed in both populations
receiver <- "New"
expressed_genes_receiver <- get_expressed_genes(receiver, seurat_obj = seuratObj, pct = 0.1, assay_oi = "RNA")
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

sender <- c("V", "P", "F1", "F2", "F3", "E", "Neuts", "Mono", "Macs", "cDC1", "cDC2", "NK", "Tcells", "Bcells")
list_expressed_genes_sender <- sender %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.1, "RNA")
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

#2. Define a gene set of interest: these are the genes in the "receiver/target" cell population tht are potentially affected by ligands expressed by interacting cells (e.g genes differentially expressed upon cell-cell interaction)
seurat_obj_receiver <- subset(seuratObj, idents = receiver)
seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["aggregate"]])

condition_oi <- "Repairing"
condition_reference <- "Control"

DE_table_receiver <- FindMarkers(object=seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct=0.1) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#3. Define a set of potential ligands: these are ligands that are expressed by the "sender/niche" cell population and bind a (putative) receptor exressed by the "receiver/target" population
ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#4. Perform NicheNet ligand activity analysis: rank the potential ligands based on the presence of their target genes in the gene set of interest (compared to the background set of genes)
ligand_activities <- predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = bckground_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities <- ligand_activities %>% arrange(-pearson) %>% mutate(rank=rank(desc(pearson)))
ligand_activities 
# the activity measures (auroc, aupr, pearson correlation coefficient) are a measure for how well a ligand can predict the observed DEgenes compared tothe background of expressed genes
# pearson correlatoin coefficient between a ligands target predictions and the observed transcriptional response was hte most informative measure to define ligand activity
# prioritise ligands inducing DEgene changes in New cells
best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
#visualise expression of top ligands by sender cell populations
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu")+ RotatedAxis()
VlnPlot(seuratObj, features=best_upstream_ligands, ncol = 5, pt.size = 0)

#5. Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
#active target gene inference to see which top-ranked ligands are predicted to have regulated the expression of which differentially expressed genes
active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset= geneset_oi, ligand_target_matrix= ligand_target_matrix, n=200) %>% bind_rows() %>% drop_na()
active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()

#visualise
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0090, 0.02))
p_ligand_target_network

#dotplot
DotPlot(seuratObj %>% subset(idents="New"), features=order_targets %>% rev(), split.by = "aggregate")+RotatedAxis()
DotPlot(seuratObj, features = order_receptors %>% rev(), cols = "RdYlBu")+ RotatedAxis()
#Vlnplot
VlnPlot(seuratObj %>% subset(idents="New"), features= c("Pcdh7", "Cebpb", "Fn1", "Mmp10", "Mmp14", "Mmp3", "Tnc", "Sod2"), ncol=4, split.by = "aggregate", pt.size = 0)
#featureplot
FeaturePlot(seuratObj, features= c("Pcdh7", "Cebpb", "Fn1", "Mmp10"), ncol=4, split.by = "aggregate", pt.size = 0.5, min.cutoff = 0)
FeaturePlot(seuratObj, features= c("Mmp14", "Mmp3", "Tnc", "Sod2"), ncol=4, split.by = "aggregate", pt.size = 0.5, min.cutoff = 0)

#receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
#visualise
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#dotplot
DotPlot(seuratObj %>% subset(idents="New"), features=order_receptors %>% rev(), split.by = "aggregate")+RotatedAxis()
DotPlot(seuratObj, features = order_receptors %>% rev(), cols = "RdYlBu")+ RotatedAxis()
#Vlnplot
VlnPlot(seuratObj %>% subset(idents="New"), features= c("Pcdh7", "Itgav", "Itgb3", "Sdc4", "Fzd1", "Il1r1", "Bdkrb2", "S1pr2"), ncol=4, split.by = "aggregate", pt.size = 0)
#featureplot
FeaturePlot(seuratObj, features= c("Pcdh7", "Itgav", "Itgb3", "Sdc4"), ncol=4, split.by = "aggregate", pt.size = 0.5, min.cutoff = 0)
FeaturePlot(seuratObj, features= c("Fzd1", "Il1r1", "Bdkrb2", "S1pr2"), ncol=4, split.by = "aggregate", pt.size = 0.5, min.cutoff = 0)

#receptors of top-ranked ligands, but after considering only bona fide liagnd-receptor interactions documented in literaure and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
#visulise
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

#6. add log fold change information of ligands from sended cells
# in some cases it is possible to check upregulation of ligands in sender cells. this can dd a useful extra lyer of information next to the ligand activites defined by nichenet because you cn assume that some of the ligands inducing DE in receiver cells, will be DE themselves in the sender cells
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "aggregate", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10) %>% reduce(full_join)
DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

#7. Summary visualisation of the NicheNet analysis
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))

# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype %in% sender), features = order_ligands_adapted %>% rev(), cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) 

figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

FeaturePlot(seuratObj, reduction = "umap", features = "Ackr3")


#############################################
library(SingleCellSignalR)

seuratObj[["integrated"]]@counts


nichenet_output <- nichenet_seuratobj_aggregate(receiver = receiver, seuratObj, condition_colname = "group", condition_oi = "Repairing", condition_reference = "Control", sender= sender, ligand_target_matrix, lr_network, weighted_networks_lr)
