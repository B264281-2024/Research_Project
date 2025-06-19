#load packages - 20240731
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(CellChat)
source("Utils.R")
library(tidyverse)
library(tsne)
library(NMF)
library(circlize)
library(ComplexHeatmap)
library(CellChat)

#Things I had to install before CellChat

BiocManager::install("devtools")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocNeighbors")

devtools::install_github("jokergoo/ComplexHeatmap")
install.packages('NMF')
devtools::install_github("jokergoo/circlize")
library(devtools)

#Installing CellChat

devtools::install_github("sqjin/CellChat")

#Setting work directory on mac
setwd("~/OneDrive - University of Edinburgh/PhD/Single cell/2023_RJA_analysis")

#Read in final data file 

data <- readRDS(file = "2023_RJA_leukocytes_combined_final.Rds")

#-------------------------------------------------------------------------------
#20240731 

#Workflow below comes from the following CellChat tutorial: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html#preprocessing-the-expression-data-for-cell-cell-communication-analysis

#PART I: Data input and processing and initialization of CellChat object
#extract cellchat input files from Seurat object 
data <- readRDS(file = "2023_RJA_leukocytes_combined_final.Rds")
Idents(data) <- data$celltype_2
DimPlot(data)
data.input <- GetAssayData(data, assay = "RNA", slot = "data")

labels <- Idents(data)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use1 <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

#set the used database in the object
cellchat@DB <- CellChatDB
cellchat@DB

#preprocessing the expression data for cell-cell communication analysis. #To infer cell state-specific communications, need to identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat, only.pos = TRUE, thresh.pc = 0.10)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)

#PART II: inference of cell-cell communication network 
#compute the communication probability and infer cellular communication network

options(future.globals.maxSize = 1000 * 1024^2)
        
cellchat <- computeCommunProb(cellchat)

#compute the communication probability and infer cellular communication network but set raw.use = FALSE 
cellchat_FALSE <- computeCommunProb(cellchat, trim = 0.1, raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat_FALSE <- filterCommunication(cellchat, min.cells = 10) #Doing FALSE didnt change the number of pairs so don't bother with this going forward.

saveRDS(cellchat, file = "20241102_cellchat.Rds")

#This is how far I got so read the file back in and continue from here...

cellchat <- readRDS(file = "20241030_cellchat.Rds")


#extract the inferred cellular communication network as a data frame 
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
write.csv(df.net, file="cell_cell_communication.csv")
df.net <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.net, file="pathways_communication.csv")

#infer cell-cell communication at a signalling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#calculate the agregated cell-cell communication netwrok 
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")

#examine signals sent from each cell group 
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Part III: Visualization of cell-cell communication network

cellchat@netP$pathways
pathways.show.all <- cellchat@netP$pathways


#Here we take input of one signaling pathway as an example. All the signaling pathways showing significant communications can be accessed by cellchat@netP$pathways.

pathways.show <- c("CCL")  
pathways.show <- c("THBS")
pathways.show <- c("SPP1")
pathways.show <- c("COMPLEMENT")
pathways.show <- c("TNF")
pathways.show <- c("BST2")

pathways.show <- c("GALECTIN")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)


pathways.show <- c("Spp1")  

selected_pair <- "Ccl9_ccr1"

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = selected_pair, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

# Plot the aggregated cell-cell communication network at the signaling pathway level

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

#Automatically save the plots of all the inferred networks for quick exploration

# Access all the signaling pathways showing significant communications

pathways.show.all <- cellchat@netP$pathways

# check the order of cell identity to set suitable vertex.receiver

levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)}

#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#Bubble plot

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1,3:4), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(5:11), remove.isolate = FALSE)


#Comparing communications on a single object

#show all the significant interactions (L-R pairs) associated with certain signaling pathways

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1,3:4), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)


#Comparing communications on a single object

#Chord diagram

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from cell type of interest
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(1,3:4), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = 2, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)


# show all the interactions received by cell type of interest
netVisual_chord_gene(cellchat, sources.use = c(1,3,4), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

#Plot the signaling gene expression distribution using violin/dot plot

#We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
plotGeneExpression(cellchat, signaling = "CXCL")

#By default, plotGeneExpression only shows the expression of signaling genes related to the inferred significant communications. USERS can show the expression of all signaling genes related to one signaling pathway by
plotGeneExpression(cellchat, signaling = "FN1", enriched.only = FALSE)

#Alternatively, USERS can extract the signaling genes related to the inferred L-R pairs or signaling pathway using extractEnrichedLR, and then plot gene expression using Seurat package.

#Part IV: Systems analysis of cell-cell communication network

#Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

pathways.show.all <- cellchat@netP$pathways

pathways.show <- c("THBS")
pathways.show <- c("")

netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(cellchat.cycling, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.repair, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.remodel, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
#We also provide another intutive way to visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
#We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

ggsave("TRIAL.png", plot = ht1, width = 10, height = 5, units = "in", dpi = 300)

install.packages("BiocManager")  # If not installed
BiocManager::install("ComplexHeatmap")  # Install ComplexHeatmap if missing

library(ComplexHeatmap)  # Load the package

png("/Users/rebeccaainslie/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Single cell/2023_RJA_analysisTRIAL.png", width = 10, height = 10, units = "in", res = 300)
draw(ht1)  # Use draw() for ComplexHeatmap objects
dev.off()

class(ht1)

file.exists("/Users/rebeccaainslie/Library/CloudStorage/OneDrive-UniversityofEdinburgh/PhD/Single cell/2023_RJA_analysisTRIAL.png")




#End of tutorial -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#PART I: Data input and processing and initilisation of CellChat object
#extract cellchat inpit files from Seurat object 

library(ggalluvial)

data <- readRDS(file = "2023_RJA_leukocytes_combined_final.Rds")
Idents(data) <- "group"

cycling <- subset(data, idents="Cycling")
repair <- subset(data, idents="24hrs")
remodel <- subset(data, idents="48hrs")

Idents(cycling) <- "celltype_2"
Idents(repair) <- "celltype_2"
Idents(remodel) <- "celltype_2"

data.input.cycling <- GetAssayData(cycling, assay = "RNA", slot = "data")
data.input.repair <- GetAssayData(repair, assay = "RNA", slot = "data")
data.input.remodel <- GetAssayData(remodel, assay = "RNA", slot = "data")

labels <- Idents(cycling)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

labels <- Idents(repair)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

labels <- Idents(remodel)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

cellchat.cycling <- createCellChat(object = data.input.cycling, meta = meta, group.by = "labels")
cellchat.repair <- createCellChat(object = data.input.repair, meta = meta, group.by = "labels")
cellchat.remodel <- createCellChat(object = data.input.remodel, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use1 <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

#set the used database in the object
cellchat.cycling@DB <- CellChatDB
cellchat.repair@DB <- CellChatDB
cellchat.remodel@DB <- CellChatDB
cellchat.cycling@DB

#PRE-PROCESSING EACH DATASET BEFORE MERGING
#preprocessing the expression data for cell-cell communication analysis. #To infer cell state-specific communications, need to identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
# subset the expression data of signaling genes for saving computation cost
cellchat.cycling <- subsetData(cellchat.cycling) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchat.cycling <- identifyOverExpressedGenes(cellchat.cycling, only.pos = TRUE, thresh.pc = 0.10)
cellchat.cycling <- identifyOverExpressedInteractions(cellchat.cycling)
# project gene expression data onto PPI network (optional)
cellchat.cycling <- projectData(cellchat.cycling, PPI.mouse)

#PART II: inference of cell-cell communication network 
#compute the communication probability and infer cellular communiction network
cellchat.cycling <- computeCommunProb(cellchat.cycling, trim = 0.25)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.cycling <- filterCommunication(cellchat.cycling, min.cells = 10)

#extract the inferred cellular communication network as a data frame 
df.net <- subsetCommunication(cellchat.cycling) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
write.csv(df.net, file="cell_cell_communication_cycling.csv")
df.net <- subsetCommunication(cellchat.cycling, slot.name = "netP")
write.csv(df.net, file="pathways_communication_cycling.csv")

#infer cell-cell communication at a signalling pathway level
cellchat.cycling <- computeCommunProbPathway(cellchat.cycling)

#calculate the agregated cell-cell communication netwrok 
cellchat.cycling <- aggregateNet(cellchat.cycling)
groupSize <- as.numeric(table(cellchat.cycling@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat.cycling@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.cycling@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")

#examine signals sent from each cell group 
mat <- cellchat.cycling@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

##PART IV: systems analysis of cell-cell communication
# Compute the network centrality scores
cellchat.cycling <- netAnalysis_computeCentrality(cellchat.cycling, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- c("CCL")
netAnalysis_signalingRole_network(cellchat.cycling, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
library(NMF)
library(ggalluvial)

#outgoing
selectK(cellchat.cycling, pattern = "outgoing")
nPatterns = 4
cellchat.cycling <- identifyCommunicationPatterns(cellchat.cycling, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat.cycling, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat.cycling, pattern = "outgoing")

#incoming
selectK(cellchat.cycling, pattern = "incoming")
nPatterns = 3
cellchat.cycling <- identifyCommunicationPatterns(cellchat.cycling, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat.cycling, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat.cycling, pattern = "incoming")

saveRDS(cellchat.remodel, file = "cellchat_remodel_leukocytes.rds")
saveRDS(cellchat.cycling, file = "cellchat_cycling_leukocytes.rds")
saveRDS(cellchat.repair, file = "cellchat_repair_leukocytes.rds")

cellchat.remodel <- readRDS(file = "cellchat_remodel_leukocytes.rds")

#=============================================================================================================================================================

cellchat.cycling <- readRDS(file = "cellchat_cycling_leukocytes.rds")
cellchat.repair <- readRDS(file = "cellchat_repair_leukocytes.rds")
cellchat.remodel <- readRDS(file = "cellchat_remodel_leukocytes.rds")

#merge datasets
object.list <- list(Cycling=cellchat.cycling.NMM, Repair=cellchat.repair.NMM, Remodel=cellchat.remodel.NMM)

cellchat2 <- mergeCellChat(object.list, add.names=names(object.list), cell.prefix = TRUE)
cellchat2

gg1 <- compareInteractions(cellchat2, show.legend = F, group = c(1:3))
gg2 <- compareInteractions(cellchat2, show.legend = F, group = c(1:3), measure = "weight")
gg1 + gg2

#Compare the number of interactions and interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat2, weight.scale = T)
netVisual_diffInteraction(cellchat2, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat2)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat2, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Differential number of interactions or interaction strength among different cell types
#We then can show the number of interactions or interaction strength between any two cell types in each dataset.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. ## Identify signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat2, idents.use = "Neutrophils", comparison = c(1, 2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat2, idents.use = "Neutrophils", comparison = c(1, 3))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat2, idents.use = "Neutrophils", comparison = c(2, 3))
patchwork::wrap_plots(plots = list(gg1,gg2, gg3))

gg1 + gg2 + gg3

#Part II: Identify the conserved and context-specific signaling pathways
#Identify signaling networks with larger (or less) difference as well as signaling groups based on their functional/structure similarity
#Identify signaling groups based on their functional similarity
future::plan("multisession", workers = 1)

cellchat2 <- computeNetSimilarityPairwise(cellchat2, type = "functional")
cellchat2 <- netEmbedding(cellchat2, type = "functional")
cellchat2 <- netClustering(cellchat2, type = "functional")
netVisual_embeddingPairwise(cellchat2, type = "functional", label.size = 3.5)

cellchat2 <- computeNetSimilarityPairwise(cellchat2, type = "structural")
cellchat2 <- netEmbedding(cellchat2, type = "structural")
cellchat2 <- netClustering(cellchat2, type = "structural")
netVisual_embeddingPairwise(cellchat2, type = "structural", label.size = 3.5)

#Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat2, type = "functional")

#Identify and visualize the conserved and context-specific signaling pathways
#Compare the overall information flow of each signaling pathway
gg1 <- rankNet(cellchat2, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2, 3))
gg2 <- rankNet(cellchat2, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(2, 3))
gg1 + gg2

#Compare outgoing (or incoming) signaling associated with each cell population
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 15)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 15)

draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 + ht2 + ht3

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 15, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 15, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i+2]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+2], width = 6, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#Identify dysfunctional signaling by comparing the communication probabities
netVisual_bubble(cellchat2, sources.use = 2:3, targets.use = c(2),  comparison = c(1, 2, 3), angle.x = 45, signaling = "CCL")
#> Comparing communications on a merged object

#Moreover, we can identify the upgulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(cellchat2, sources.use = 1:3, targets.use = c(1),  comparison = c(1, 2, 3), signaling= "CCL", max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat2, sources.use = 1:3, targets.use = c(2:3),  comparison = c(1, 2, 3), signaling= "CCL", max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
gg3 <- netVisual_bubble(cellchat2, sources.use = 1:3, targets.use = c(3),  comparison = c(1, 2, 3), signaling= "CCL", max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
gg1 + gg2 + gg3

#Identify dysfunctional signaling by using differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Repair"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat2 <- identifyOverExpressedGenes(cellchat2, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat2, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat2, net = net, datasets = "Repair",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat2, net = net, datasets = "Repair",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat2)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat2)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat2, pairLR.use = pairLR.use.up, sources.use = 1:3, targets.use = c(2), comparison = c(1, 2, 3), signaling = "CCL",  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat2, pairLR.use = pairLR.use.down, sources.use = 1:3, targets.use = c(2), comparison = c(1, 2, 3), signaling = "CCL", angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

#Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
pathways.show <- c("CCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1,2,3), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))

# Chord diagram
par(mfrow = c(1,2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

#Part V: Compare the signaling gene expression distribution between different datasets
cellchat2@meta$datasets = factor(cellchat2@meta$datasets, levels = c("Cycling", "Repair", "Remodel")) # set factor level
plotGeneExpression(cellchat2, signaling = "CCL", split.by = "datasets", colors.ggplot = T)

saveRDS(cellchat2, file = "cellchat_leukocytes_merged.rds")

#================================================================================================================================================================================
#Pheeb's script which I have modified for my own purpose.

cellchat.cycling <- readRDS(file = "cellchat_cycling_leukocytes.rds")
cellchat.repair <- readRDS(file = "cellchat_repair_leukocytes.rds")
cellchat.remodel <- readRDS(file = "cellchat_remodel_leukocytes.rds")

cellchat.cycling.NMM <- readRDS(file = "cellchat_cycling_NMM.rds")
cellchat.repair.NMM <- readRDS(file = "cellchat_repair_NMM.rds")
cellchat.remodel.NMM <- readRDS(file = "cellchat_remodel_NMM.rds")


#Identify overexpressed Ligand-Receptor genes
cellchat.cycling <- identifyOverExpressedGenes(cellchat.cycling)
cellchat.cycling <- identifyOverExpressedInteractions(cellchat.cycling)
cellchat.repair <- identifyOverExpressedGenes(cellchat.repair)
cellchat.repair <- identifyOverExpressedInteractions(cellchat.repair)
cellchat.remodel <- identifyOverExpressedGenes(cellchat.remodel)
cellchat.remodel <- identifyOverExpressedInteractions(cellchat.remodel)

#THIS PART OF THE SCRIPT DOESNT WORK -> either figure out how to fix it out figure out how to hide certain cell types in the plot
cellchat.cycling.NMM <- identifyOverExpressedGenes(cellchat.cycling.NMM)
cellchat.cycling.NMM <- identifyOverExpressedInteractions(cellchat.cycling.NMM)
cellchat.repair.NMM <- identifyOverExpressedGenes(cellchat.repair.NMM)
cellchat.repair.NMM <- identifyOverExpressedInteractions(cellchat.repair.NMM)
cellchat.remodel.NMM <- identifyOverExpressedGenes(cellchat.remodel.NMM)
cellchat.remodel.NMM <- identifyOverExpressedInteractions(cellchat.remodel.NMM)

#Compute communication probability an infer cellular communication network
cellchat.cycling <- computeCommunProb(cellchat.cycling, trim = 0.1)
cellchat.repair <- computeCommunProb(cellchat.repair, trim = 0.1)
cellchat.remodel <- computeCommunProb(cellchat.remodel, trim = 0.1)

#Infer cell-cell communication at a signalling pathway level, stored in the new data slots 'net' and 'netP'
cellchat.cycling <- computeCommunProbPathway(cellchat.cycling)
cellchat.repair <- computeCommunProbPathway(cellchat.repair)
cellchat.remodel <- computeCommunProbPathway(cellchat.remodel)

cellchat.cycling <- aggregateNet(cellchat.cycling)
cellchat.repair <- aggregateNet(cellchat.repair)
cellchat.remodel <- aggregateNet(cellchat.remodel)
groupSize1 <- as.numeric(table(cellchat.cycling@idents))
groupSize2 <- as.numeric(table(cellchat.cycling@idents))
groupSize3 <- as.numeric(table(cellchat.remodel@idents))
par(mfrow = c(1,2), xpd = TRUE)

netVisual_circle(cellchat.cycling@net$count, vertex.weight = groupSize1, weight.scale = T, label.edge = F, title.name = "Number of interactions: Cycling")
netVisual_circle(cellchat.repair@net$count, vertex.weight = groupSize2, weight.scale = T, label.edge = F, title.name = "Number of interactions: Repair")
netVisual_circle(cellchat.remodel@net$count, vertex.weight = groupSize3, weight.scale = T, label.edge = F, title.name = "Number of interactions: Remodel")

netVisual_circle(cellchat.cycling@net$weight, vertex.weight = groupSize1, weight.scale = T, label.edge = F, title.name = "Weight/Strength of interactions: Cycling")
netVisual_circle(cellchat.repair@net$weight, vertex.weight = groupSize2, weight.scale = T, label.edge = F, title.name = "Weight/Strength of interactions: Repair")
netVisual_circle(cellchat.remodel@net$weight, vertex.weight = groupSize3, weight.scale = T, label.edge = F, title.name = "Weight/Strength of interactions: Remodel")

#Next visualise the number and weight of communications by each cell type to one another
mat1 <- cellchat.cycling@net$weight
mat2 <- cellchat.repair@net$weight
mat3 <- cellchat.remodel@net$weight
par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat1)) {
  mat1.2 <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
  mat1.2[i, ] <- mat1[i, ]
  netVisual_circle(mat1.2, vertex.weight = groupSize1, weight.scale = T, edge.weight.max = max(mat1), title.name = rownames(mat1)[i])
}

par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat2)) {
  mat2.2 <- matrix(0, nrow = nrow(mat2), ncol = ncol(mat2), dimnames = dimnames(mat2))
  mat2.2[i, ] <- mat2[i, ]
  netVisual_circle(mat2.2, vertex.weight = groupSize2, weight.scale = T, edge.weight.max = max(mat2), title.name = rownames(mat2)[i])
}

par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat3)) {
  mat3.2 <- matrix(0, nrow = nrow(mat3), ncol = ncol(mat3), dimnames = dimnames(mat3))
  mat3.2[i, ] <- mat3[i, ]
  netVisual_circle(mat3.2, vertex.weight = groupSize3, weight.scale = T, edge.weight.max = max(mat3), title.name = rownames(mat3)[i])
}

#Show all the signalling pathways showing significant communication, returns via the console as a list
cellchat.cycling@netP$pathways
cellchat.repair@netP$pathways
cellchat.remodel@netP$pathways

cellchat.cycling <- netAnalysis_computeCentrality(cellchat.cycling, slot.name = "netP")
cellchat.repair <- netAnalysis_computeCentrality(cellchat.repair, slot.name = "netP")
cellchat.remodel <- netAnalysis_computeCentrality(cellchat.remodel, slot.name = "netP")

#Create heatmaps to show which are dominant senders and receivers for each significant signalling pathway
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.cycling, pattern = "outgoing", title = "Cycling - outgoing", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Reds")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.cycling, pattern = "incoming", title = "Cycling - incoming", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Greens")
ht1 
ht2

mt1 <- netAnalysis_signalingRole_heatmap(cellchat.repair, pattern = "outgoing", title = "Repair - outgoing", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Reds")
mt2 <- netAnalysis_signalingRole_heatmap(cellchat.repair, pattern = "incoming", title = "Repair - incoming", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Greens")
mt1
mt2

tt1 <- netAnalysis_signalingRole_heatmap(cellchat.remodel, pattern = "outgoing", title = "Remodelling - outgoing", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Reds")
tt2 <- netAnalysis_signalingRole_heatmap(cellchat.remodel, pattern = "incoming", title = "Remodelling - incoming", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Greens")
tt1
tt2

#Next ID global communication patterns, show significant dropoff for bpth cophenetic and silhouette for both cyling and T48 to determine k
library(ggalluvial)
library(NMF)
selectK(cellchat.cycling, pattern = "incoming")
selectK(cellchat.repair, pattern = "incoming")
selectK(cellchat.remodel, pattern = "incoming")
selectK(cellchat.cycling, pattern = "outgoing")
selectK(cellchat.repair, pattern = "outgoing")
selectK(cellchat.remodel, pattern = "outgoing")

#Make a heatmap to visualise the patterns via cell type and signalling pathways
cellchat.cycling <- identifyOverExpressedGenes(cellchat.cycling)
cellchat.repair <- identifyOverExpressedGenes(cellchat.repair)
cellchat.remodel <- identifyOverExpressedGenes(cellchat.remodel)

cellchat.cycling <- identifyCommunicationPatterns(cellchat.cycling, pattern = "incoming", k = 4)
cellchat.repair <- identifyCommunicationPatterns(cellchat.repair, pattern = "incoming", k = 4)
cellchat.remodel <- identifyCommunicationPatterns(cellchat.remodel, pattern = "incoming", k = 4)

cellchat.cycling <- identifyCommunicationPatterns(cellchat.cycling, pattern = "outgoing", k = 4)
cellchat.repair <- identifyCommunicationPatterns(cellchat.repair, pattern = "outgoing", k = 4)
cellchat.remodel <- identifyCommunicationPatterns(cellchat.remodel, pattern = "outgoing", k = 4)

#Next make river plots to show better
netAnalysis_river(cellchat.cycling, pattern = "incoming")
netAnalysis_river(cellchat.repair, pattern = "incoming")
netAnalysis_river(cellchat.remodel, pattern = "incoming")

netAnalysis_river(cellchat.cycling, pattern = "outgoing")
netAnalysis_river(cellchat.repair, pattern = "outgoing")
netAnalysis_river(cellchat.remodel, pattern = "outgoing")

#Can use dot plots to visualise relative contributions to each pattern
netAnalysis_dot(cellchat.cycling, pattern = "incoming")
netAnalysis_dot(cellchat.repair, pattern = "incoming")
netAnalysis_dot(cellchat.remodel, pattern = "incoming")

netAnalysis_dot(cellchat.cycling, pattern = "outgoing")
netAnalysis_dot(cellchat.repair, pattern = "outgoing")
netAnalysis_dot(cellchat.remodel, pattern = "outgoing")

#Next want to identify signalling groups based on either structural or functional similarity
reticulate::py_install(packages ='umap-learn')

#Can't get this bit to work

cellchat.cycling <- computeNetSimilarity(cellchat.cycling, type = "functional")
cellchat.cycling <- netEmbedding(cellchat.cycling, type = "functional")
cellchat.cycling <- netClustering(cellchat.cycling, type = "functional")
netVisual_embedding(cellchat.cycling, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat.cycling, type = "functional", nCol = 2)

cellchat.repair <- computeNetSimilarity(cellchat.repair, type = "functional")
cellchat.repair <- netEmbedding(cellchat.repair, type = "functional")
cellchat.repair <- netClustering(cellchat.repair, type = "functional")
netVisual_embedding(cellchat.repair, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat.repair, type = "functional", nCol = 2)

cellchat.remodel <- computeNetSimilarity(cellchat.remodel, type = "functional")
cellchat.remodel <- netEmbedding(cellchat.remodel, type = "functional")
cellchat.remodel <- netClustering(cellchat.remodel, type = "functional")
netVisual_embedding(cellchat.remodel, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat.remodel, type = "functional", nCol = 2)

#Chord plot per pathway, combine similar clusters together
pathways.show1 <- c("THBS")
par(mfrow = c(1,1))
netVisual_aggregate(cellchat.cycling, signaling = pathways.show1, layout = "circle")

pathways.show2 <- c("CXCL")
par(mfrow = c(1,1))
netVisual_aggregate(cellchat.repair, signaling = pathways.show2, layout = "circle")

pathways.show1 <- c("THBS")
par(mfrow = c(1,1))
netVisual_aggregate(cellchat.remodel, signaling = pathways.show1, layout = "circle")

#Can analyse the contributions of each ligand-receptor pair to the overall signalling pathway
netAnalysis_contribution(cellchat.cycling, signaling = pathways.show1)
netAnalysis_contribution(cellchat.repair, signaling = pathways.show1)
netAnalysis_contribution(cellchat.remodel, signaling = pathways.show1)

#Can also visualise cell-cell communication mediated by a single L-R pair
PairLR <- extractEnrichedLR(cellchat.cycling, signaling = pathways.show1, geneLR.return = FALSE)
LR.show <- PairLR[1,3]
vertex.receiver = seq(1,4)
netVisual_individual(cellchat.cycling, signaling = pathways.show1, pairLR.use = LR.show, vertex.receiver = vertex.receiver)

PairLR <- extractEnrichedLR(cellchat.remodel, signaling = pathways.show1, geneLR.return = FALSE)
LR.show <- PairLR[1,3]
vertex.receiver = seq(1,4)
netVisual_individual(cellchat.remodel, signaling = pathways.show3, pairLR.use = LR.show, vertex.receiver = vertex.receiver)

PairLR <- extractEnrichedLR(cellchat.repair, signaling = pathways.show2, geneLR.return = FALSE)
LR.show <- PairLR[1,3]
vertex.receiver = seq(1,4)
netVisual_individual(cellchat.repair, signaling = pathways.show2, pairLR.use = LR.show, vertex.receiver = vertex.receiver)

cycling_lr <- netVisual_individual(cellchat.cycling, signaling = pathways.show1, pairLR.use = LR.show, layout = "circle")
repair_lr <- netVisual_individual(cellchat.repair, signaling = pathways.show1, pairLR.use = LR.show, layout = "circle")
remodel_lr <- netVisual_individual(cellchat.remodel, signaling = pathways.show1, pairLR.use = LR.show, layout = "circle")


#Create Violin plots to establish expression levels of the unique identified LR genes for each celltype cluster 
cv1 <- VlnPlot(cycling, features = c("Ccl6", "Ccr1"), pt.size = 0, assay = "RNA") + scale_y_continuous(limits = c(0,5))
tv1 <- VlnPlot(remodel, features = c("Plxnd1"), pt.size = 0, assay = "RNA") + scale_y_continuous(limits = c(0,5))
cv1 + tv1

cp1 <- DimPlot(cycling, reduction = "umap", label = FALSE, repel = TRUE, group.by = "orig.ident", label.size = 5, pt.size = 0.5)
cp2 <- DimPlot(cycling, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5, pt.size = 0.5)  
cp2 + cp1

tp1 <- DimPlot(remodel, reduction = "umap", label = FALSE, repel = TRUE, group.by = "orig.ident", label.size = 5, pt.size = 0.5)
tp2 <- DimPlot(remodel, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5, pt.size = 0.5)
tp2 + tp1

tp3 <- FeaturePlot(remodel, features = "Cxcl3", pt.size = 1)
tp3
tp2

#Analysis of single cell RNA seq data using Seurat
#using nichnetR fr ligand-receptor anlaysis

#data format: table of read or UMI counts, one per column per individual cell and one row per gene
BiocManager::install ("Seurat")
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)
library(pheatmap)
theme_set(theme_cowplot())

#read in seurat objects that you wish to integrate generated using the seurat_training_script.R
data1 <- readRDS(file="seurat_object1")
data2 <- readRDS(file="seurat_object2")

#add metadata to describe which cells this dataset represents
data1@meta.data$group <- rep("Celltype1", times=length(data1))
data2@meta.data$group <- rep("Celltype2", times=length(data2))

#merge two objects
data.merged <- merge(data1, y=data2, add.cell.ids = c("celltype1", "celltype2"), project = "control cells1/cells2 merge")
head(colnames(data.merged))
table(data.merged$orig.ident)
data.merged

#pre-process new merged seurat object
data.merged <- ScaleData(data.merged)
data.merged <- FindVariableFeatures(object=data.merged)
data.merged <- RunPCA(data.merged, npcs=30, features = VariableFeatures(object = data.merged))
data.merged <- RunUMAP(data.merged, reduction= "pca", dims=1:20)
data.merged <- FindNeighbors(data.merged, dims = 1:20)
data.merged <- FindClusters(data.merged, resolution=0.1)
head(Idents(data.merged), 5)

DimPlot(data, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 6, group.by = "seurat cluster") #assign cell types and change idents and oder as required

saveRDS(data, file = "celltpe1_celltype2_merged_Seurat_Object.rds")

#====================================================================================================================================================================================================
#Using NicheNetR to perform ligand-receptor analysis

#I started from here with cellchat2 object (all leukocytes from all experimental conditions). 

# install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")

library(nichenetr)
library(tidyverse)

#check data
seuratObj <- readRDS(file = "2023_RJA_leukocytes_combined_final.Rds")
DimPlot(seuratObj, reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 6, group.by = "group")
seuratObj@meta.data %>% head()

#visualise which cell types are present in the data
seuratObj@meta.data$celltype_2 %>% table()
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

#====================================================================================================================================================================================================
##PERFORM THE NICHENET ANALYSIS
#apply NicheNet to predict which ligands expressed by all cells in the dataset are most likely to have induced the differential expression in CELLS OF INTEREST after condition change (ie during repair)
#NicheNet analysis pipeline
#1. Define a "sender/niche" cell population and a "receiver/target" cell population present in your expression data and determine which genes are expressed in both populations

Idents(object = seuratObj) <- "celltype_2"

receiver <- c("Neutrophils")
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.1, assay_oi = "RNA")
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

sender <- c("Macrophages", "Monocytes")
list_expressed_genes_sender <- sender %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.1, "RNA")
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

expressed_genes_receiver
expressed_genes_sender

#2. Define a gene set of interest: these are the genes in the "receiver/target" cell population tht are potentially affected by ligands expressed by interacting cells (e.g genes differentially expressed upon cell-cell interaction)
seurat_obj_receiver <- subset(seuratObj, idents = receiver)

seurat_obj_receiver <- SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["group"]]) #CANT GET THIS LINE TO WORK

Idents(object = seuratObj) <- "group"

condition_oi <- "24hrs"
condition_reference <- "48hrs"

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

head(geneset_oi)
geneset_oi
#3-define a set of potential ligands which are expressed by sender population and bind to a putative receptor on receiver population 
ligands = lr_network %>% pull(from) %>% unique()
ligands
receptors = lr_network %>% pull(to) %>% unique()
receptors
expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_ligands
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands
#4perform nicheet ligand activity analysis: rank potential ligands based on the presence of their target genes in the gene set of interest 9compaed with background genes)
#auroc, aupr, pearson: measures for how well a ligand can predict observed differentilly expressed genes compared to background - allows prioritisation of ligands 
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

seuratObj@meta.data %>% head()

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(seuratObj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

#5 infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis 
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
#visualise: 

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

#dotplot
DotPlot(seuratObj, features=order_targets %>% rev(), split.by = "group", cols = "RdYlBu")+RotatedAxis()

#Vlnplot
Idents(object = seuratObj)
Idents(object = seuratObj) <- "celltype_2"
VlnPlot(seuratObj, features= c("Spp1", "Vegfa", "Csf1", "Lgals3", "Ccl2"), ncol=5, pt.size = 0) 
mono_mac <- subset(x=seuratObj, idents="Monocytes/macrophages")
fibro_peri<- subset( x=seuratObj, idents =  c("Fibroblasts", "Perivascular cells"))
fibro<- subset(x=seuratObj, idents = "Fibroblasts")
Idents(object = fibro)<- "group"
Idents(object = fibro_peri) <- "group"
Idents(object = mono_mac) <- "group"
VlnPlot(mono_mac, features= c("Spp1", "Vegfa", "Csf1", "Ccl2", "Lgals3"), ncol=5, split.by = "group", pt.size = 0) 
VlnPlot(fibro_peri, features= c("Col5a3", "Bmp7", "Pgf", "Ctgf", "Mmp13"), ncol=5, split.by = "group", pt.size = 0) 
VlnPlot(fibro_peri, features= c("Tgfb"), ncol=5, split.by = "group", pt.size = 0) 
#"Timp1", "Mmp3", "Col1a1", "Mmp13"
#featureplot
FeaturePlot(seuratObj, features= c("Spp1", "Vegfa", "Csf1", "Ccl2", "Lgals3"), ncol=5, split.by = "group", pt.size = 0.5, min.cutoff = 0)


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

#visualise data 
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#receptors of top-ranked ligands but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
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

#visualising expression of receptors on mac/mono across groups 

Idents(object = seuratObj) <- "celltype_2"

DotPlot(seuratObj, features= order_receptors %>% rev(), split.by = "group", cols = "RdYlBu")+RotatedAxis() + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) 
VlnPlot(seuratObj, features= c("Cd44", "Tnfrsf1a", "Tnfrsf1b", "Notch1", "Il1r2", "Ilrap", "Fpr1", "Fpr2", "Sorl1", "Icam1", "Alox5"), ncol=4, split.by = "group", pt.size = 0) 

genes2<- c("Cd44", "Tnfrsf1a", "Tnfrsf1b", "Notch1", "Il1r2", "Il1rap")
genes3 <- c("Fpr1", "Fpr2", "Sorl1", "Icam1", "Alox5")
StackedVlnPlot(seuratObj, features = genes2)
StackedVlnPlot(seuratObj, features = genes3)

#d log fold chnage information from sender cells. check upreglation of ligands in sender cells. To test if ligands causing DE in rceiving cells are all themseleves have DE. 
# DE analysis for each sender cell type
# this uses a new nichenetr function - reinstall nichenetr if necessary!
DE_table_all = Idents(seuratObj) %>% levels() %>% intersect(sender) %>% lapply(get_lfc_celltype, seurat_obj = seuratObj, condition_colname = "group", condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = 0.10, celltype_col = "celltype_2") %>% reduce(full_join) # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
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

# change colors a bit to make them more stand out
p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
p_ligand_lfc

#summary of data visualtion of whole analysis 
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson
# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
rotated_dotplot = DotPlot(seuratObj %>% subset(celltype_2 %in% sender), features = order_ligands_adapted %>% rev(), cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) 
rotated_dotplot
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

#20241118 RJA Make a CellChat object for each of the experimental conditions but only containing neuts, monos and macs

data <- readRDS(file = "2023_RJA_leukocytes_combined_final.Rds")
Idents(data) <- "group"

cycling <- subset(data, idents="Cycling")
repair <- subset(data, idents="24hrs")
remodel <- subset(data, idents="48hrs")

Idents(cycling) <- "celltype_2"
Idents(repair) <- "celltype_2"
Idents(remodel) <- "celltype_2"

cycling_NMM <- subset(cycling, idents=c("Neutrophils", "Monocytes", "Macrophages"))
repair_NMM <- subset(repair, idents=c("Neutrophils", "Monocytes", "Macrophages"))
remodel_NMM <- subset(remodel, idents=c("Neutrophils", "Monocytes", "Macrophages"))

Idents(cycling_NMM)
Idents(repair_NMM)
Idents(remodel_NMM)

data.input.cycling.NMM <- GetAssayData(cycling_NMM, assay = "RNA", slot = "data")
data.input.repair.NMM <- GetAssayData(repair_NMM, assay = "RNA", slot = "data")
data.input.remodel.NMM <- GetAssayData(remodel_NMM, assay = "RNA", slot = "data")

labels <- Idents(remodel_NMM)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

cellchat.cycling.NMM <- createCellChat(object = data.input.cycling.NMM, meta = meta, group.by = "labels")
cellchat.repair.NMM <- createCellChat(object = data.input.repair.NMM, meta = meta, group.by = "labels")
cellchat.remodel.NMM <- createCellChat(object = data.input.remodel.NMM, meta = meta, group.by = "labels")

saveRDS(cellchat.remodel.NMM, file = "cellchat_remodel_NMM.rds")
saveRDS(cellchat.cycling, file = "cellchat_cycling_NMM.rds")
saveRDS(cellchat.repair, file = "cellchat_repair_NMM.rds")


#=================================================================================================================

#3/2/24: Plots for manuscript

plotGeneExpression(cellchat, signaling = "SPP1", enriched.only = FALSE)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.cycling <- netAnalysis_computeCentrality(cellchat.cycling, slot.name = "netP")
cellchat.repair <- netAnalysis_computeCentrality(cellchat.repair, slot.name = "netP")
cellchat.remodel <- netAnalysis_computeCentrality(cellchat.remodel, slot.name = "netP")

netAnalysis_signalingRole_network(cellchat, signaling = "THBS", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.cycling, signaling = "THBS", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.repair, signaling = "THBS", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.remodel, signaling = "THBS", width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(cellchat, signaling = "SPP1", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.cycling, signaling = "SPP1", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.repair, signaling = "SPP1", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.remodel, signaling = "SPP1", width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(cellchat, signaling = "MIF", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.cycling, signaling = "MIF", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.repair, signaling = "MIF", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.remodel, signaling = "MIF", width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(cellchat, signaling = "CCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.cycling, signaling = "CCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.repair, signaling = "CCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.remodel, signaling = "CCL", width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(cellchat, signaling = "CXCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.cycling, signaling = "CXCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.repair, signaling = "CXCL", width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat.remodel, signaling = "CXCL", width = 8, height = 2.5, font.size = 10)

#11/2/25: Run some of cell chat just with NMM.

Neuts_Monos_Macs <- readRDS(file = "Neuts_Monos_Macs.Rds")

data.input <- GetAssayData(Neuts_Monos_Macs, assay = "RNA", slot = "data")

labels <- Idents(Neuts_Monos_Macs)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use1 <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

#set the used database in the object
cellchat@DB <- CellChatDB
cellchat@DB

#preprocessing the expression data for cell-cell communication analysis. #To infer cell state-specific communications, need to identify over-expressed ligands or receptors in one cell group and then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat, only.pos = TRUE, thresh.pc = 0.10)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)

#PART II: inference of cell-cell communication network 
#compute the communication probability and infer cellular communication network

options(future.globals.maxSize = 1000 * 1024^2)

cellchat <- computeCommunProb(cellchat)

#compute the communication probability and infer cellular communication network but set raw.use = FALSE 
cellchat_FALSE <- computeCommunProb(cellchat, trim = 0.1, raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat_FALSE <- filterCommunication(cellchat, min.cells = 10) #Doing FALSE didnt change the number of pairs so don't bother with this going forward.

saveRDS(cellchat, file = "20250211_cellchat_NMM.Rds")

#This is how far I got so read the file back in and continue from here...

cellchat_NMM <- readRDS(file = "20250211_cellchat_NMM.Rds")

#extract the inferred cellular communication network as a data frame 
df.net <- subsetCommunication(cellchat_NMM) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
write.csv(df.net, file="cell_cell_communication.csv")
df.net <- subsetCommunication(cellchat_NMM, slot.name = "netP")
write.csv(df.net, file="pathways_communication.csv")

#infer cell-cell communication at a signalling pathway level
cellchat_NMM <- computeCommunProbPathway(cellchat_NMM)

#calculate the agregated cell-cell communication netwrok 
cellchat_NMM <- aggregateNet(cellchat_NMM)
groupSize <- as.numeric(table(cellchat_NMM@idents))
par(mfrow = c(1,1), xpd=TRUE)
netVisual_circle(cellchat_NMM@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_NMM@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weight/strength")

#examine signals sent from each cell group 
mat <- cellchat_NMM@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Part III: Visualization of cell-cell communication network

cellchat_NMM@netP$pathways
pathways.show.all <- cellchat_NMM@netP$pathways


#Here we take input of one signaling pathway as an example. All the signaling pathways showing significant communications can be accessed by cellchat_NMM@netP$pathways.

pathways.show <- c("CCL")  
pathways.show <- c("THBS")
pathways.show <- c("SPP1")
pathways.show <- c("COMPLEMENT")
pathways.show <- c("TNF")
pathways.show <- c("BST2")

pathways.show <- c("APP")

# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 

vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_NMM, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_NMM, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat_NMM, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat_NMM, signaling = pathways.show, color.heatmap = "Reds")

# Do heatmap based on a single object

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat_NMM@idents)
netVisual_chord_cell(cellchat_NMM, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

# Plot the aggregated cell-cell communication network at the signaling pathway level

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair

netAnalysis_contribution(cellchat_NMM, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat_NMM, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[4,] # show one ligand-receptor pair

# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat_NMM, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat_NMM, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

#Chord diagram
netVisual_individual(cellchat_NMM, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

#Automatically save the plots of all the inferred networks for quick exploration

# Access all the signaling pathways showing significant communications

pathways.show.all <- cellchat_NMM@netP$pathways

# check the order of cell identity to set suitable vertex.receiver

levels(cellchat_NMM@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat_NMM, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat_NMM, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)}

#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#Bubble plot

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat_NMM, sources.use = 1, targets.use = c(1,3:4), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 2, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 3, targets.use = c(1,3:4), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 4, targets.use = c(1,3:4), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 5, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 6, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 7, targets.use = c(5:11), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 1, targets.use = c(1), remove.isolate = FALSE)

#Comparing communications on a single object

#show all the significant interactions (L-R pairs) associated with certain signaling pathways

netVisual_bubble(cellchat_NMM, sources.use = 1, targets.use = c(1:3), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 2, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 3, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 5, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 6, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 7, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
netVisual_bubble(cellchat_NMM, sources.use = 8, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)

#Comparing communications on a single object

#Chord diagram

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from cell type of interest
netVisual_chord_gene(cellchat_NMM, sources.use = 1, targets.use = c(1:3), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_NMM, sources.use = 2, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_NMM, sources.use = 3, targets.use = c(1,3:4), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat_NMM, sources.use = 4, targets.use = c(1,3:4), lab.cex = 0.5,legend.pos.y = 30)

# show all the interactions received by cell type of interest
netVisual_chord_gene(cellchat_NMM, sources.use = c(1,2,3), targets.use = 8, legend.pos.x = 15)

# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat_NMM, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat_NMM, sources.use = c(1,2,3,4), targets.use = c(5:11), slot.name = "netP", legend.pos.x = 10)

#Plot the signaling gene expression distribution using violin/dot plot

#We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
plotGeneExpression(cellchat_NMM, signaling = "CXCL")

#By default, plotGeneExpression only shows the expression of signaling genes related to the inferred significant communications. USERS can show the expression of all signaling genes related to one signaling pathway by
plotGeneExpression(cellchat_NMM, signaling = "FN1", enriched.only = FALSE)

#Alternatively, USERS can extract the signaling genes related to the inferred L-R pairs or signaling pathway using extractEnrichedLR, and then plot gene expression using Seurat package.

#Part IV: Systems analysis of cell-cell communication network

#Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat_NMM <- netAnalysis_computeCentrality(cellchat_NMM, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

pathways.show.all <- cellchat_NMM@netP$pathways

pathways.show <- c("THBS")
pathways.show <- c("")

netAnalysis_signalingRole_network(cellchat_NMM, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(cellchat_NMM.cycling, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat_NMM.repair, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat_NMM.remodel, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)



#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
#We also provide another intutive way to visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat_NMM)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat_NMM, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
#We can also answer the question on which signals contributing most to outgoing or incoming signaling of certain cell groups.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_NMM, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_NMM, pattern = "incoming")
ht1 + ht2


#20250211: Just using NMM for some of the cell chat plots

#20241204 Generating figures for aim 1 manuscript:

Neuts_Monos_Macs <- readRDS(file = "Neuts_Monos_Macs.Rds")

Idents(Neuts_Monos_Macs) <- Neuts_Monos_Macs$group
Neuts_Monos_Macs_Cyc <- subset(x=Neuts_Monos_Macs, idents = c("Cycling")) 
Neuts_Monos_Macs_Rep <- subset(x=Neuts_Monos_Macs, idents = c("24hrs")) 
Neuts_Monos_Macs_Res <- subset(x=Neuts_Monos_Macs, idents = c("48hrs")) 

Neuts_Monos_Macs_Cyc <- RunUMAP(Neuts_Monos_Macs_Cyc, reduction="pca", dims=1:20)
Neuts_Monos_Macs_Rep <- RunUMAP(Neuts_Monos_Macs_Rep, reduction="pca", dims=1:20)
Neuts_Monos_Macs_Res <- RunUMAP(Neuts_Monos_Macs_Res, reduction="pca", dims=1:20)

Neuts_Monos_Macs_Cyc <- FindNeighbors(Neuts_Monos_Macs_Cyc, dims = 1:20)
Neuts_Monos_Macs_Rep <- FindNeighbors(Neuts_Monos_Macs_Rep, dims = 1:20)
Neuts_Monos_Macs_Res <- FindNeighbors(Neuts_Monos_Macs_Res, dims = 1:20)

DimPlot(Neuts_Monos_Macs_Cyc, reduction = "umap", pt.size = 0.5, label = FALSE, label.size = 6)
DimPlot(Neuts_Monos_Macs_Rep, reduction = "umap", pt.size = 0.5, label = FALSE, label.size = 6)
DimPlot(Neuts_Monos_Macs_Res, reduction = "umap", pt.size = 0.5, label = FALSE, label.size = 6)

saveRDS(Neuts_Monos_Macs_Cyc, file = "Cycling_NMM.Rds")
saveRDS(Neuts_Monos_Macs_Rep, file = "Repair_NMM.Rds")
saveRDS(Neuts_Monos_Macs_Res, file = "Remodel_NMM.Rds")

Idents(Neuts_Monos_Macs) <- Neuts_Monos_Macs$celltype_2

#PART I: Data input and processing and initialization of CellChat object
#extract cellchat input files from Seurat object 
Neuts_Monos_Macs_Cyc <- readRDS(file = "Cycling_NMM.Rds")
Neuts_Monos_Macs_Rep <- readRDS(file = "Repair_NMM.Rds")
Neuts_Monos_Macs_Res <- readRDS(file = "Remodel_NMM.Rds")

Idents(Neuts_Monos_Macs_Cyc) <- Neuts_Monos_Macs_Cyc$celltype_2
Idents(Neuts_Monos_Macs_Rep) <- Neuts_Monos_Macs_Rep$celltype_2
Idents(Neuts_Monos_Macs_Res) <- Neuts_Monos_Macs_Res$celltype_2

data.input.cycling.NMM <- GetAssayData(Neuts_Monos_Macs_Cyc, assay = "RNA", slot = "data")
data.input.repair.NMM <- GetAssayData(Neuts_Monos_Macs_Rep, assay = "RNA", slot = "data")
data.input.remodel.NMM <- GetAssayData(Neuts_Monos_Macs_Res, assay = "RNA", slot = "data")

labels <- Idents(Neuts_Monos_Macs_Res)
meta = data.frame(group = labels, row.names = names(labels))
colnames(meta) <- "labels"

cellchat.cycling.NMM <- createCellChat(object = data.input.cycling.NMM, meta = meta, group.by = "labels")
cellchat.repair.NMM <- createCellChat(object = data.input.repair.NMM, meta = meta, group.by = "labels")
cellchat.remodel.NMM <- createCellChat(object = data.input.remodel.NMM, meta = meta, group.by = "labels")

#THIS PART OF THE SCRIPT DOESNT WORK -> either figure out how to fix it out figure out how to hide certain cell types in the plot
# Managed to get cycling to work but not repair or remodelling?
cellchat.cycling.NMM <- identifyOverExpressedGenes(cellchat.cycling.NMM)
cellchat.cycling.NMM <- identifyOverExpressedInteractions(cellchat.cycling.NMM)
cellchat.repair.NMM <- identifyOverExpressedGenes(cellchat.repair.NMM)
cellchat.repair.NMM <- identifyOverExpressedInteractions(cellchat.repair.NMM)
cellchat.remodel.NMM <- identifyOverExpressedGenes(cellchat.remodel.NMM)
cellchat.remodel.NMM <- identifyOverExpressedInteractions(cellchat.remodel.NMM)

#Compute communication probability an infer cellular communication network
cellchat.cycling.NMM <- computeCommunProb(cellchat.cycling.NMM, trim = 0.1)
cellchat.repair <- computeCommunProb(cellchat.repair, trim = 0.1)
cellchat.remodel <- computeCommunProb(cellchat.remodel, trim = 0.1)

#Infer cell-cell communication at a signalling pathway level, stored in the new data slots 'net' and 'netP'
cellchat.cycling.NMM <- computeCommunProbPathway(cellchat.cycling.NMM)
cellchat.repair <- computeCommunProbPathway(cellchat.repair)
cellchat.remodel <- computeCommunProbPathway(cellchat.remodel)

cellchat.cycling.NMM <- aggregateNet(cellchat.cycling.NMM)
cellchat.repair <- aggregateNet(cellchat.repair)
cellchat.remodel <- aggregateNet(cellchat.remodel)
groupSize1 <- as.numeric(table(cellchat.cycling.NMM@idents))
groupSize2 <- as.numeric(table(cellchat.repair@idents))
groupSize3 <- as.numeric(table(cellchat.remodel@idents))
par(mfrow = c(1,2), xpd = TRUE)

netVisual_circle(cellchat.cycling.NMM@net$count, vertex.weight = groupSize1, weight.scale = T, label.edge = F, title.name = "Number of interactions: Cycling")
netVisual_circle(cellchat.repair@net$count, vertex.weight = groupSize2, weight.scale = T, label.edge = F, title.name = "Number of interactions: Repair")
netVisual_circle(cellchat.remodel@net$count, vertex.weight = groupSize3, weight.scale = T, label.edge = F, title.name = "Number of interactions: Remodel")

netVisual_circle(cellchat.cycling.NMM@net$weight, vertex.weight = groupSize1, weight.scale = T, label.edge = F, title.name = "Weight/Strength of interactions: Cycling")
netVisual_circle(cellchat.repair@net$weight, vertex.weight = groupSize2, weight.scale = T, label.edge = F, title.name = "Weight/Strength of interactions: Repair")
netVisual_circle(cellchat.remodel@net$weight, vertex.weight = groupSize3, weight.scale = T, label.edge = F, title.name = "Weight/Strength of interactions: Remodel")

#Next visualise the number and weight of communications by each cell type to one another
mat1 <- cellchat.cycling.NMM@net$weight
mat2 <- cellchat.repair@net$weight
mat3 <- cellchat.remodel@net$weight
par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat1)) {
  mat1.2 <- matrix(0, nrow = nrow(mat1), ncol = ncol(mat1), dimnames = dimnames(mat1))
  mat1.2[i, ] <- mat1[i, ]
  netVisual_circle(mat1.2, vertex.weight = groupSize1, weight.scale = T, edge.weight.max = max(mat1), title.name = rownames(mat1)[i])
}

par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat2)) {
  mat2.2 <- matrix(0, nrow = nrow(mat2), ncol = ncol(mat2), dimnames = dimnames(mat2))
  mat2.2[i, ] <- mat2[i, ]
  netVisual_circle(mat2.2, vertex.weight = groupSize2, weight.scale = T, edge.weight.max = max(mat2), title.name = rownames(mat2)[i])
}

par(mfrow = c(3,4), xpd = TRUE)
for (i in 1:nrow(mat3)) {
  mat3.2 <- matrix(0, nrow = nrow(mat3), ncol = ncol(mat3), dimnames = dimnames(mat3))
  mat3.2[i, ] <- mat3[i, ]
  netVisual_circle(mat3.2, vertex.weight = groupSize3, weight.scale = T, edge.weight.max = max(mat3), title.name = rownames(mat3)[i])
}

#Show all the signalling pathways showing significant communication, returns via the console as a list
cellchat.cycling.NMM@netP$pathways
cellchat.repair@netP$pathways
cellchat.remodel@netP$pathways

cellchat.cycling.NMM <- netAnalysis_computeCentrality(cellchat.cycling.NMM, slot.name = "netP")
cellchat.repair <- netAnalysis_computeCentrality(cellchat.repair, slot.name = "netP")
cellchat.remodel <- netAnalysis_computeCentrality(cellchat.remodel, slot.name = "netP")

#Create heatmaps to show which are dominant senders and receivers for each significant signalling pathway
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.cycling.NMM, pattern = "outgoing", title = "Cycling - outgoing", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Reds")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.cycling.NMM, pattern = "incoming", title = "Cycling - incoming", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Greens")
ht1 
ht2

mt1 <- netAnalysis_signalingRole_heatmap(cellchat.repair, pattern = "outgoing", title = "Repair - outgoing", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Reds")
mt2 <- netAnalysis_signalingRole_heatmap(cellchat.repair, pattern = "incoming", title = "Repair - incoming", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Greens")
mt1
mt2

tt1 <- netAnalysis_signalingRole_heatmap(cellchat.remodel, pattern = "outgoing", title = "Remodelling - outgoing", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Reds")
tt2 <- netAnalysis_signalingRole_heatmap(cellchat.remodel, pattern = "incoming", title = "Remodelling - incoming", font.size = 11, font.size.title = 13, height = 18, color.heatmap = "Greens")
tt1
tt2

#merge datasets
object.list <- list(Cycling=cellchat.cycling.NMM, Repair=cellchat.repair.NMM, Remodel=cellchat.remodel.NMM)

cellchat_NMM <- mergeCellChat(object.list, add.names=names(object.list), cell.prefix = TRUE)
cellchat_NMM

gg1 <- compareInteractions(cellchat_NMM, show.legend = F, group = c(1:3))
gg2 <- compareInteractions(cellchat_NMM, show.legend = F, group = c(1:3), measure = "weight")
gg1 + gg2

#Compare the number of interactions and interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_NMM, weight.scale = T)
netVisual_diffInteraction(cellchat_NMM, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat_NMM)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_NMM, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

#To better control the node size and edge weights of the inferred networks across different datasets, we compute the maximum number of cells per cell group and the maximum number of interactions (or interaction weights) across all datasets.
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Differential number of interactions or interaction strength among different cell types
#We then can show the number of interactions or interaction strength between any two cell types in each dataset.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Compare the major sources and targets in 2D space

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

#Furthermore, we can identify the specific signaling changes of Inflam.DC and cDC1 between NL and LS. ## Identify signaling changes associated with one cell group
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_NMM, idents.use = "Neutrophils", comparison = c(1, 2))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_NMM, idents.use = "Neutrophils", comparison = c(1, 3))
gg3 <- netAnalysis_signalingChanges_scatter(cellchat_NMM, idents.use = "Neutrophils", comparison = c(2, 3))
patchwork::wrap_plots(plots = list(gg1,gg2, gg3))

gg1 + gg2 + gg3




saveRDS(cellchat.cycling.NMM, file = "cellchat_repair_NMM.Rds")
saveRDS(cellchat.rep.NMM, file = "cellchat_repair_NMM.Rds")
saveRDS(cellchat.res.NMM, file = "cellchat_remodel_NMM.Rds")