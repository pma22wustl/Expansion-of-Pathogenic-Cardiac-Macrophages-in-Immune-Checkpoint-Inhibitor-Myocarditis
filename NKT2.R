library(dplyr)
library(Seurat)
library(patchwork)
seuratobj <- readRDS("NKT_rm11_10_8.rds")
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",split.by="orig.ident",label=TRUE,pt.size = 1.5)

DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)
#Plotting multiple genes in the same FeaturePlot
marker_gene_list <- list(c("Cd8a", "Cd4", "Tyrobp"))
object <- AddModuleScore(object, features = marker_gene_list, name = "Cluster1_score1")
head(object$Cluster1_score11)
FeaturePlot(object = object, features = "Cluster1_score11",pt.size = 1.5,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE)
+
  theme(plot.title = element_text(size=10))



# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time choose a approxipriate dim

DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)

seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seuratobj), 10)
plot1 <- VariableFeaturePlot(seuratobj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:20, nfeatures = 5)
DimPlot(seuratobj, reduction = "pca")
seuratobj <- JackStraw(seuratobj, num.replicate = 100)
seuratobj <- ScoreJackStraw(seuratobj, dims = 1:20)
JackStrawPlot(seuratobj, dims = 1:15)
ElbowPlot(seuratobj)

seuratobj <- FindNeighbors(seuratobj, dims = 1:10)
seuratobj <- FindClusters(seuratobj, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0))

seuratobj <- RunUMAP(seuratobj,dims = 1:10)

#for umap
DefaultAssay(seuratobj) <- "integrated"
Idents(seuratobj) <- seuratobj$integrated_snn_res.0.3
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",split.by="orig.ident",label=TRUE,pt.size = 1.5)
saveRDS(seuratobj, file = "NKT_final_highres.rds")

library(ggplot2)
ggsave("res1.0 dim10 split.png")
DimPlot(seuratobj, reduction = "umap",label=TRUE,pt.size = 1.5)
ggsave(".png")
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,split.by = "group")

#for UMAPs
DefaultAssay(seuratobj) <- "integrated"
Idents(seuratobj) <- seuratobj$integrated_snn_res.1
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",split.by="group",label=TRUE,pt.size = 1.5)
library(ggplot2)
ggsave("res1.0 split.png")

# remove cluster
seuratobj <- subset(seuratobj, idents =  c("0","1","2","3","4","5","6","7","9"))
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 2.0)

saveRDS(seuratobj, file = "NKT_rm11_10_8.rds")
library(ggplot2)
ggsave("res1.0 dim35 split.png")
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 1.5)
ggsave("res1.png",width = 15,height = 10)



#switch assay to RNA Expression analysis

DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)

NKT.markers <- FindAllMarkers(seuratobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(NKT.markers,"NKT_markerhighres.csv")

#vlnplot
VlnPlot(seuratobj, features = c("Tnfsf4","Tnfrsf4"),split.by="group",split.plot=FALSE,pt.size = 0.0,ncol = 2)
library(ggplot2)
ggsave(filename ="genes_Tnfsf4vlnplot4_NKT.png",path="Whole object/",width = 6,height = 5)

#Featureplot 
FeaturePlot(seuratobj, features = c("Tnfsf4","Tnfrsf4","Cd69","Tnfrsf9"),pt.size = 0.5,split.by="group",label = FALSE,
            blend = FALSE,cols = c("lightgrey","red"))

ggsave("OX40.png",path = "Whole object/",width = 6,height=8)
#Featureplot two genes in one plot
FeaturePlot(seuratobj, features = c("Cxcl9", "Ccr2"),pt.size = 0.5,label = T,
            blend = TRUE )
library(ggplot2)
ggsave("Ccr2 Cxcl9featureplot.png")
library(font)
install.packages("font")