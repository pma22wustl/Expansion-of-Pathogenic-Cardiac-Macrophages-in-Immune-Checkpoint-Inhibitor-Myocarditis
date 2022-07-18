#cell-cell interatcion
if(!require(devtools)) install.packages("devtools");
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE,force = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap",force = TRUE)
devtools::install_github("jokergoo/circlize")
library(circlize)
#install.packages('NMF')
library(NMF)
remotes::install_github("sqjin/CellChat",force = TRUE)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(patchwork)
library(tidyverse)
#Part I: Data input & processing and initialization of CellChat object

#extract matrix with cell type info from two seurat project
seuratobj <- readRDS("Named_T_only.rds")
seuratobj2 <- readRDS("Onlymac_final.rds")
seuratobj.combined <- merge(seuratobj, y = seuratobj2, add.cell.ids = c("T", "myeloid"), project = "combined")
seuratobj<-seuratobj.combined

DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE) #split.by = "group"
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)
data.input  <- seuratobj@assays$RNA@data
identity = data.frame(group =seuratobj@active.ident,row.names = names(seuratobj@active.ident)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat <- createCellChat(data.input)
cellchat
summary(cellchat)
#用思维导图可视化数据结构,可省略
#install.packages("mindr")
library(mindr)
(out <- capture.output(str(cellchat)))
out2 <- paste(out, collapse="\n")
mm(gsub("\\.\\.@","# ",gsub("\\.\\. ","#",out2)),type ="text")

#把metadata信息加到CellChat对象中
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

saveRDS(cellchat,"cellchat-t&m.rds")

#将不同的 CellChat 对象合并在一起
cellchat_T <- readRDS("cellchat-t.rds")
cellchat_myeloid <- readRDS("cellchat-m.rds")
object.list <- list(T = cellchat_T, myeloid = cellchat_myeloid)
cellchat<-cellchat_myeloid
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
summary(cellchat)
cellchat

#用CellChatDB.human,CellChatDB.mouse来导入配受体数据库
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

#optional
(out3 <- capture.output(str(CellChatDB)))
out4 <- paste(out3, collapse="\n")
mm(gsub("\\$","# ",gsub("\\.\\. ","#",out4)),type ="text")

#选特定的信息描述细胞间的相互作用，即从特定的侧面来刻画细胞互作
unique(CellChatDB$interaction$annotation)
# use Secreted Signaling for cell-cell communication analysis

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use 

#在一个细胞组中识别过表达的配体或受体，然后将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。
#如果配体或受体过表达，则识别过表达配体和受体之间的相互作用

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)
# do parallel
future::plan("multiprocess", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)  

#Part II: Inference of cell-cell communication network
#通过为每个相互作用分配一个概率值并进行置换检验来推断生物意义上的细胞-细胞通信
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Extract the inferred cellular communication network as a data frame
#access the inferred cell-cell communications of interest
df.net <- subsetCommunication(cellchat)#slot.name = "netP", to access the inferred communications at the level of signaling pathways
write.csv(df.net,"cell-cell communications mac&T_cells.csv")
#get cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5
df.net <- subsetCommunication(cellchat, sources.use = c(2,3), targets.use = c(4))
#cell-cell communications mediated by signaling WNT and TGFb.
df.net <- subsetCommunication(cellchat, signaling = c("IFN-II")) 


#inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot 'net' and 'netP', respectively.
#Infer the cell-cell communication at a signaling pathway level
#Calculate the aggregated cell-cell communication network
cellchat <- computeCommunProbPathway(cellchat)
#sources.use= and targets.use=: calculate the aggregated network among a subset of cell groups
cellchat <- aggregateNet(cellchat)

# visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
#creating a simple multi-paneled plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
ggsave("communication.pdf",width=15,height = 5)
#examine the signaling sent from each cell group. control the parameter edge.weight.max to compare edge weights between differet networks.
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

#Part III: Visualization of cell-cell communication network
#Hierarchy plot: define vertex.receiver, the target cell groups,
#solid and open circles represent source and target, respectively
#Chord diagram: netVisual_chord_cell(cell groups), 
# netVisual_chord_gene(pathway, ligand or receptor


netVisual_chord_cell(2)

#查看结果
cellchat@netP$pathways
head(cellchat@LR$LRsig)
#查看细胞族群
levels(cellchat@idents) 
#visualize the signal pathway
pathways.show <- c("IFN-II") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,6) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
# Chord diagram
group.cellType <- c(rep("T", 3),rep("myeloid", 3),"T") # grouping cell clusters into T and myeloid cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

#contribution of each ligand-receptor pair to the overall signaling pathway
netAnalysis_contribution(cellchat, signaling = pathways.show)
#extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CD45 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CD45[3,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,6) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
# Chord diagram
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,6)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

saveRDS(cellchat,"cellchat myeloid&T combined after calculation.rds")
cellchat<-readRDS("cellchat myeloid&T combined after calculation.rds")

#Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
#show (L-R pairs) from some cell groups to other cell groups using netVisual_bubble
netVisual_bubble(cellchat, sources.use =c(1:4), targets.use = c(6:10), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("TGFb","PTPRM"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MHC-II","PTPRM"))
netVisual_bubble(cellchat, sources.use = c(5:9), targets.use = c(1:4), pairLR.use = pairLR.use, remove.isolate = TRUE)

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:10), slot.name = "netP", legend.pos.x = 10)
netVisual_chord_gene(cellchat, sources.use = c(5:10), targets.use =c(1,2,3,4) , slot.name = "netP", legend.pos.x = 10)

#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "CCR2")#enriched.only = FALSE show all genes
plotGeneExpression(cellchat, signaling = "CD45", )
# extract the signaling genes related to the inferred L-R pairs or signaling pathway
gene_expression<-extractEnrichedLR(cellchat,signaling = pathways.show.all)
gene_expression

#Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(4,2))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all, width = 8, height = 2.5, font.size = 10)

#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTPRM", "MHC-II"))
gg1 + gg2
gg1

 # which signals contributing most to outgoing or incoming signaling of certain cell groups.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

#Identify and visualize outgoing communication pattern of secreting cells
library(NMF)
library(ggalluvial)
#run selectK to infer the number of patterns.
selectK(cellchat, pattern = "outgoing")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

#Identify and visualize incoming communication pattern of target cells
selectK(cellchat, pattern = "incoming")
nPatterns = 4 #Cophenetic values begin to drop when the number of incoming patterns is 4.
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
install.packages("Miniconda")
install.packages("Python")
installed.packages('umap-learn')
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

#Identify signaling groups based on structure similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

saveRDS(cellchat,"cellchat myeloid&T combined after calculation.rds")
sessionInfo()





#write.table(seuratobj@assays[["RNA"]]@counts, file='Gene_Count_per_T_Cell.tsv', quote=FALSE, sep='\t', col.names = TRUE)
#T_data <- read.table(file = "Gene_Count_per_T_Cell.tsv", header = TRUE,as.is = TRUE, row.names = 1)
T_data1<-as.dataframe(t(T_data))
celltype<- as.data.frame(seuratobj@active.ident)
T_cell_matrix<-data.frame(cbind(T,celltype))

colnames(T_cell_matrix)[33695]="cell_type"
T_cell_matrix %>% remove_rownames() %>% column_to_rownames(var="cell_type")
column_to_rownames(T_cell_matrix, "cell_type")

rownames(T_cell_matrix) <- T_cell_matrix[,-1]
T_cell_matrix<-T_cell_matrix[,-1]
write.csv(T_cell_matrix, "T_cell_matrix with cell type.csv")

library(iTALK)
T_data<-read.csv("T_cell_matrix with cell type.csv",row.names = 1)#avoid the auto generated column with number of rows
myeloid_data<-read.csv("myeloid_cell_matrix with cell type.csv",row.names = 1)
mergedata<-rbind(T_data,myeloid_data)
write.csv(mergedata,"matrix for T-myeloid crosstalk.csv")

data<-mergedata[,-1]
colnames(data)[1]="symbol"
colnames(dat)[2]="avg_log2FC"
colnames(data)[-1]="cell_type"


data<-T_cell_matrix
## highly expressed ligand-receptor pairs
# find top 50 percent highly expressed genes
highly_exprs_genes<-rawParse(T_data1,top_genes=50,stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_col<-structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),names=unique(data$cell_type))
par(mfrow=c(1,2))
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
  res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
  #plot by ligand category
  #overall network plot
  NetView(res_cat,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
  #top 20 ligand-receptor pairs
  LRPlot(res_cat[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res_cat$cell_from_mean_exprs[1:20],link.arr.width=res_cat$cell_to_mean_exprs[1:20])
  title(comm_type)
  res<-rbind(res,res_cat)
}
res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),][1:20,]
NetView(res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=res$cell_from_mean_exprs[1:20],link.arr.width=res$cell_to_mean_exprs[1:20])