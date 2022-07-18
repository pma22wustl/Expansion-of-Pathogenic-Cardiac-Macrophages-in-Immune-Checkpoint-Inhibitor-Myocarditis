library(dplyr)
library(Seurat)
library(patchwork)
seuratobj <- readRDS("myeloid_rm89_reclustered.rds")
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 1.5)
ggsave("myeloid res0.4.png")
Idents(seuratobj)
levels(seuratobj)
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)
# install monocle
BiocManager::install("monocle")
library(monocle)
packageVersion("BiocManager")

library(devtools)
devtools::install_github("cole-trapnell-lab/DDRTree", ref="simple-ppt-like")
devtools::install_github("cole-trapnell-lab/L1-graph")
install.packages("reticulate")
library(reticulate)
py_install('umap-learn', pip = T, pip_ignore_installed = T) # Ensure the latest version of UMAP is installed
py_install("louvain")

devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle)
packageVersion("monocle")

#monocle构建CDS需要3个矩阵：expr.matrix、pd、fd
# 将Seurat中的对象转换为monocle识别的对象
#cds <- importCDS(GetAssayData(seuratobj)) it didn't work
#Load Seurat object
seurat_object <- readRDS('myeloid_rm89_reclustered.rds')

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

#查看phenodata、featuredata
head(pData(monocle_cds))
head(fData(monocle_cds))


#预处理
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
#筛选基因,这里可以根据自己的需要筛选特定的基因
disp_table <- dispersionTable(monocle_cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
#用DDRtree 进行降维分析
monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  method = 'DDRTree')
#计算psudotime值
monocle_cds <- orderCells(monocle_cds)
head(pData(monocle_cds))
save(monocle_cds,file="myeloid_trajectoray.Rdata")

plot_cell_trajectory(monocle_cds,cell_size = 1)

plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",cell_size = 2)
plot_cell_trajectory(monocle_cds, color_by = "State",cell_size = 2)
plot_cell_trajectory(monocle_cds, color_by = "seurat_clusters",cell_size = 1)
plot_cell_trajectory(monocle_cds, color_by = "group",cell_size = 1.0,group_by="group")
plot_cell_trajectory(monocle_cds, color_by = "integrated_snn_res.0.4",cell_size =2.0)
plot_cells

plot_cells(cds=monocle_cds,color_cells_by="Pseudotime",show_trajectory_graph=TRUE
           )

library(ggplot2)
ggsave("myeloid0.4psedotime1.0 group.png")
