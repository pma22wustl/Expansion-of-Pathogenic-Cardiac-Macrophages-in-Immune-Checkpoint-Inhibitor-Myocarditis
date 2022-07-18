library(dplyr)
library(Seurat)
library(patchwork)
seuratobj <- readRDS("myeloid_rm89_reclustered.rds")
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,split.by = "group")

seuratobj[["percent.rb"]] <- PercentageFeatureSet(seuratobj2, pattern = "^Rp[sl]")

seuratobj$log10GenesPerUMI <- log10(seuratobj$nFeature_RNA)/log10(seuratobj2$nCount_RNA)
# remove cluster0 
seuratobj <- subset(seuratobj, ident= c("0","1","2","3","4","5","7"))

DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,split.by = "group")
saveRDS(seuratobj,"myeloid_rm6_reclustered.rds")
meyloid<-readRDS("Myeloid_rename_rm6.rds")
DimPlot(meyloid, dims=c(1,2), reduction = "umap",label=FALSE,pt.size = 1.5)
Mac<-subset(meyloid,idents = c("Cd163 resident Mac","Cxcl9 Gbp2b Mac","Nlrp3 Mac"))
DimPlot(Mac, dims=c(1,2), reduction = "umap",label=FALSE,pt.size = 1.5)
saveRDS(Mac,"Onlymac_final.rds")
Mon_Mac<-Mac<-subset(meyloid,idents = c("Cd163 resident Mac","Cxcl9 Gbp2b Mac","Nlrp3 Mac","Classic Mon"))
DimPlot(Mon_Mac, dims=c(1,2), reduction = "umap",label=FALSE,pt.size = 1.5)
saveRDS(Mon_Mac,"Mon&Mac_final.rds")
#Rename each cluster 
DefaultAssay(seuratobj) <- "integrated"
Idents(seuratobj) <- seuratobj$celltype
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,split.by = "group")
seuratobj<- RenameIdents(seuratobj, `0` = "Cd163 resident Mac", `1` = "Cxcl9 Gbp2b Mac", `2` = "Nlrp3 Mac", 
                         `3` = "Cd163 resident Mac",`4` = "Cd163 resident Mac",  `5` = "DC", `7` ="Classic Mon" )
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=FALSE,pt.size = 1.5)
seuratobj$celltype<-seuratobj@active.ident
saveRDS(seuratobj,"Myeloid_rename_rm6.rds")
seuratobj<-readRDS("Myeloid_rename_rm6.rds")
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=FALSE,pt.size = 1.5)
library(ggplot2)
ggsave("res0.4 new total with name_splitbyorig.ident.png",path = "Fianl analysis_1_19_22/",width = 10,height = 5)
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",split.by="group",label=FALSE,pt.size = 1.5)
ggsave("res0.4 new split with name.png",path = "Fianl analysis_1_19_22/")

#cluster cell population in each group
library(scales)
table(Idents(seuratobj))
table(Idents(seuratobj),seuratobj$group)
prop.table(table(Idents(seuratobj)))
cell.prop<-as.data.frame(prop.table(table(Idents(seuratobj), seuratobj$group)))
colnames(cell.prop)<-c("cluster","group","proportion")
ggplot(cell.prop,aes(group,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  scale_y_continuous(labels = percent)+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  
  guides(fill=guide_legend(title=NULL))
ggsave("percent.png",width=4.5,height=4.5,path = "Fianl analysis_1_19_22/")

# Proplot 
source('PropPlot.R')
install.packages("ggstatsplot")
install.packages("gginnards")
library(ggstatsplot)
table(seuratobj$celltype)
PropPlot(object = seuratobj, groupBy = 'celltype')  
 ggsave("group.composition.png",width=4.5,height=4.5,path = "Fianl analysis_1_19_22/")

# cluster composition
table(Idents(seuratobj),seuratobj$group)
cell.count<- as.data.frame(table(Idents(seuratobj),seuratobj$group))
colnames(cell.count)<-c("cluster","group","count")
write.csv(cell.count,"cell.count.csv")

DF_Count2<-cell.count %>% group_by(cluster)%>%
  mutate(freq = count / sum(count)*100)

#DF_Count2 <- seuratobj@meta.data %>%group_by(integrated_snn_res.0.4,group) %>%
# summarise(n = n()) %>%
#mutate(freq = n / sum(n)*100)

DF_Count2$group <- factor(DF_Count2$group, levels = c("Ctrl", "ICI"))

png(file = "Myeloid0.4_Cell Composition Percent.png", 
    width = 15*100,
    height = 10*100,
    res = 100)
ggplot(DF_Count2, aes(x =cluster, y = freq, fill = group))+
  geom_col() +
  labs(x="Groups",y="Percentage",fill="Clusters")+
  #geom_text(aes(label = paste(round(freq, 2),"%")),position = position_stack(vjust = 0.5))+
  theme_classic()
dev.off()

#switch assay to RNA Expression analysis

DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)

#plot with original data
#vlnplot
VlnPlot(seuratobj, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb","log10GenesPerUMI"),split.plot=FALSE,pt.size = 0.01,ncol = 2)
library(ggplot2)
ggsave(filename ="genes vlnplot_with_features.png",width=10,height=10,path = "Fianl analysis_1_19_22/")

#Featureplot 
FeaturePlot(seuratobj, features = c("Cd163","Fcgr4"),pt.size = 1.5,label = FALSE,
            blend = TRUE)
ggsave("Fcgr4Cd163.png",width = 15,height = 10)
#Featureplot two genes in one plot
FeaturePlot(seuratobj, features = c("Stat1"),pt.size = 1.0,label = T,
            blend = FALSE )
library(ggplot2)
ggsave("Ccr2 Cxcl9featureplot.png")
library(font)
install.packages("font")
object<-seuratobj
#Plotting multiple genes in the same FeaturePlot
marker_gene_list <- list(c("Irf5"))
object <- AddModuleScore(object, features = marker_gene_list, name = "Cluster1_score1")
head(object$Cluster1_score11)
FeaturePlot(object = object, features = "Cluster1_score11",pt.size = 2.0,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE,split.by = "group")+
  labs(title="Irf5")+
  theme(plot.title = element_text(size=10))
ggsave("Irf5_2.png",width=10, height=5,path = "Fianl analysis_1_19_22/")

marker_gene_list2 <- list(c("H2-T22","H2-Eb1","H2-Ab1","Cd74","Ctss","H2-DMb1","Psme1","Psme2","Fcgr1","B2m","Psmb9","Psap"))
object <- AddModuleScore(object, features = marker_gene_list2, name = "Cluster2_score1")
head(object$Cluster2_score11)
FeaturePlot(object = object, features = "Cluster2_score11",pt.size = 2.0,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE)+
  labs(title="Antigen processing and presention ")+
  theme(plot.title = element_text(size=10))
ggsave("Antigen_processing_2.png",width=5, height=5,path = "Fianl analysis_1_19_22/")
marker_gene_list3 <- list(c("H2-Eb1","H2-Ab1","Ccl8","Ccl5","Ccl4","Stat1","Irgm1","Gbp2b","Gbp2","Gbp7","Gbp4","Cd40","Gbp8"))
object <- AddModuleScore(object, features = marker_gene_list3, name = "Cluster3_score1")
head(object$Cluster3_score11)
FeaturePlot(object = object, features = "Cluster3_score11",pt.size = 2.0,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE)+
  labs(title="Response to interferon-gamma")+
  theme(plot.title = element_text(size=10))
ggsave("Response to interferon-gamma_2.png",width=5, height=5,path = "Fianl analysis_1_19_22/")
marker_gene_list4 <- list(c("Slamf8","Cd9","Mif","Lgals3","Aif1","Il1b","Cxcl9","Cxcl10","Ccl5","Ccl8","Cd74"))
object <- AddModuleScore(object, features = marker_gene_list4, name = "Cluster4_score1")
head(object$Cluster4_score11)
FeaturePlot(object = object, features = "Cluster4_score11",pt.size = 2.0,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE)+
  labs(title="Myeloid leukocyte migration")+
  theme(plot.title = element_text(size=10))
ggsave("Myeloid leukocyte migration_2.png",width=5,height=5,path = "Fianl analysis_1_19_22/")
library(ggplot2)
ggsave("z score DC featureplot.png")

marker_gene_list5 <- list(c("Ly6c2","Sell","Chil3","Ace","Hp"))
object <- AddModuleScore(object, features = marker_gene_list5, name = "Cluster5_score1")
head(object$Cluster5_score11)
FeaturePlot(object = object, features = "Cluster5_score11",pt.size = 1.5,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE)+
  labs(title="Ly6c2,Sell,Chil3,Ace,Hp")+
  theme(plot.title = element_text(size=10))

library(ggplot2)
ggsave("MONOCYTE_markergenes.png",path = "Fianl analysis_1_19_22/")


# Dotplot
Idents(seuratobj) <- factor(Idents(seuratobj))
markers.to.plot <- c("Cbr2", "Gbp2b","Nlrp3","Nabp1","Plac8","Clec9a" )
DotPlot(seuratobj, features = rev(markers.to.plot), cols = c("grey","red"), dot.scale = 10 
) + RotatedAxis()

# optimized marker for dotplot
Idents(seuratobj) <- factor(Idents(seuratobj))
markers.to.plot <- c("Ccr2","Cbr2","Lyve1","Cd163","Cxcl9", "Gbp2b","Cxcl10","Ccl4","Nlrp3","Flt3","Ccr7","Clec9a","Chil3","Ly6c2","Plac8")
DotPlot(seuratobj, features = rev(markers.to.plot), cols = c("white","red"), dot.scale = 6,cluster.idents = TRUE 
) + RotatedAxis()
ggsave("genesdotplot2.pdf",path = "Fianl analysis_1_19_22/")

#Ccr2 density 
BiocManager::install("Nebulosa")
library(Nebulosa)
plot_density(seuratobj,"Mta3")
library(ggplot2)
ggsave("density_Mta3.png",width=5, height=5,path = "Fianl analysis_1_19_22/")

# Find marker gene
Idents(seuratobj)<-seuratobj$group
ICI.markers <- FindMarkers(seuratobj, ident.1 = "ICI", min.pct = 0.25)
top10 <-ICI.markers %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seuratobj, features = rownames(top10)) + NoLegend()
write.csv(ICI.markers,"ICI_3_27.markers.csv") 
ggsave("heatmap.png",path = "Fianl analysis_1_19_22/")

# Fealture plot of top10 genes

marker_gene_list5 <- rownames(top10)
object <- AddModuleScore(object, features = marker_gene_list5, name = "Cluster5_score1")
head(object$Cluster5_score11)
FeaturePlot(object = object, features = "Cluster5_score11",pt.size = 1.5,cols=c("darkblue","cyan","springgreen","yellow","orange","red"),label = FALSE,split.by = "group")+
  labs(title="top10")+
  theme(plot.title = element_text(size=10))

ggsave("zscore_top10.png",path = "Fianl analysis_1_19_22/")

# Find all marker genes
myeloid.markers <- FindAllMarkers(seuratobj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myeloid.markers,"myeloid_marker CtrlvsICInew.csv")
# Heatmap 
myeloid.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- myeloid.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(seuratobj, features = top10$gene) + NoLegend()
ggsave("Heatmap 5.png",width = 15,height = 10)
#select genes to do heatmap
for_heatmap <- myeloid.markers%>% 
  filter(gene %in% c("Cbr2","Folr2","Cd163","Lyve1","Ccl24","Fcgr4","Cxcl9", "Ccl8","Gbp2b","Cxcl10","Cxcl2","Nlrp3","Bcl2a1a","Cd14","Ccl4","cd24a","Flt3","Serpinb6b","Bh1he40","Ccr7","Ly6c2","Ace","Hp","Chil3","Sell")) 
DoHeatmap(seuratobj, features = for_heatmap$gene)
ggsave("Heatmap selected gene.png")


#switch assay to integrated
DefaultAssay(seuratobj) <- "integrated"
Idents(seuratobj) <- seuratobj$integrated_snn_res.0.4
DimPlot(seuratobj, dims=c(1,2), reduction = "umap",label=TRUE,pt.size = 0.9)
library(ggplot2)
ggsave("res 0.4 com.png")

#switch assay to RNA
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)
#plot with original data
VlnPlot(seuratobj, features = c("Cxcl9"),split.by = "group",pt.size = 0.01)
FeaturePlot(seuratobj, features = c("Cxcl9"),split.by = "group",pt.size = 0.5,
            cols =  c("grey", "red") )
#combined violinplot
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)
plots <- VlnPlot(seuratobj, features = c("Cxcl9", "Ccl8","Gbp2b","Cxcl10","Fcgr4","Gbp2"), split.by = "group",
                 pt.size = 0)
CombinePlots(plots = plots, ncol = 1)
FeaturePlot(seuratobj, features = c("Cxcl9", "Ccl8","Gbp2b","Cxcl10","Fcgr4","Gbp2"), split.by = "group", max.cutoff = 3, 
            cols = c("grey", "red"))
# Find marker gene
ICI.markers <- FindMarkers(seuratobj, ident.1 = "ICI", min.pct = 0.25)
head(ICI.markers, n = 20)
write.csv(ICI.markers,"ICI_3_27.markers.csv") 

# Find all marker genes
myeloid.markers <- FindAllMarkers(seuratobj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(myeloid.markers,"myeloid_marker CtrlvsICInew.csv")
# Heatmap 
myeloid.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- myeloid.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(seuratobj, features = top10$gene) + NoLegend()
ggsave("Heatmap 5.png",width = 15,height = 10)
#select genes to do heatmap
for_heatmap <- myeloid.markers%>% 
  filter(gene %in% c("Cbr2","Folr2","Cd163","Lyve1","Ccl24","Fcgr4","Cxcl9", "Ccl8","Gbp2b","Cxcl10","Cxcl2","Nlrp3","Bcl2a1a","Cd14","Ccl4","cd24a","Flt3","Serpinb6b","Bh1he40","Ccr7","Ly6c2","Ace","Hp","Chil3","Sell")) 
DoHeatmap(seuratobj, features = for_heatmap$gene)
ggsave("Heatmap selected gene.png")

# Exported data for palantir analysis
# Save the normalized SCT matrix
Idents(seuratobj)<-seuratobj$celltype
DefaultAssay(seuratobj) <- "RNA"
seuratobj <- NormalizeData(seuratobj)
seuratobj <- ScaleData(seuratobj)

write.csv(as.matrix(seuratobj[["RNA"]]@scale.data),
          file = "F:/2019 Lavine lab/Single-cell analysis Pan ICI/NEW/myeloid_final_palantir_RNA.txt", quote = FALSE)
getwd()
# Save the meta data
write.csv(seuratobj@meta.data, file = "F:/2019 Lavine lab/Single-cell analysis Pan ICI/NEW/myeloid_final_palantir_meta.csv", quote = TRUE)

# After palantir analysis waterfall pseudotime 
library(ggplot2)
pseud<-read.csv("./Fianl analysis_1_19_22/palantir_re_analysis/palantir_meta_data_final.csv")
install.packages("waterfalls")
library(waterfalls)
colnames(pseud)
ggplot(pseud,aes(x=pseudotime,y=entropy,colour=ClusterName))+ geom_point(size=1.5)+theme_bw()+theme_classic()
ggsave("waterfall_myeloid.pdf",width = 6,height = 6,path = "Fianl analysis_1_19_22/palantir_re_analysis/")
#只选一组
s<-subset(pseud,pseud$ClusterName=="Nlrp3 Mac")
ggplot(s,aes(x=pseudotime,y=entropy,colour=ClusterName))+ geom_point(size=1.5)+theme_bw()+theme_classic()
ggsave("Nlrp3_waterfall.pdf",width = 6,height = 6,path = "Fianl analysis_1_19_22/palantir_re_analysis/")      

p<-ggplot(pseud,aes(x=pseudotime,y=entropy))+ geom_point(aes(color=factor(ClusterName)),size=1.5)+theme_bw()+theme_classic()
p+scale_color_manual(values = c("Cd163 resident Mac"="#F08080","Cxcl9 Gbp2b Mac"="#698B23","DC"="#63B8FF","Nlrp3 Mac"="#66CDAA","Classic Mon"="#E066FF"))
ggsave("waterfall_myeloid_all.pdf",width = 6,height = 3,path = "Fianl analysis_1_19_22/palantir_re_analysis/")

p+scale_color_manual(values = c("Cd163 resident Mac"="grey","Cxcl9 Gbp2b Mac"="grey","DC"="grey","Nlrp3 Mac"="#66CDAA","Classic Mon"="#E066FF"))
ggsave("waterfall_myeloid_Nlrp3.pdf",width = 6,height = 3,path = "Fianl analysis_1_19_22/palantir_re_analysis/")

# pseudotime information added to seuratobj
seuratobj@meta.data$annotation<-seuratobj@active.ident
seuratobj@meta.data$Pseudotime<-pseud$pseudotime
seuratobj@meta.data$Entropy<-pseud$entropy
seuratobj@meta.data$TCD163<-pseud$Cd163.resident.Mac
seuratobj@meta.data$TCXCL9<-pseud$Cxcl9.Gbp2b.Mac
seuratobj@meta.data$CXCL9<-seuratobj@meta.data$Cxcl9
seuratobj@meta.data<- seuratobj@meta.data[,-grep("Cxcl9",colnames(seuratobj@meta.data))]
seuratobj@meta.data$TDC<-pseud$DC
saveRDS(seuratobj,"Rename_myeloidrm6withPseudotime.rds")
seuratobj<-readRDS("Rename_myeloidrm6withPseudotime.rds")

#Pseudotime to UMAP 
library(Seurat)
library(ggplot2)
p<-FeaturePlot(seuratobj, feature=c("Pseudotime","Entropy"),pt.size = 0.5)
p<-FeaturePlot(seuratobj, feature=c("TDC"),pt.size = 0.5)
p+scale_colour_viridis_c(option = "plasma")
ggsave("TCXCL9.pdf",width=5,height=5,path = "Fianl analysis_1_19_22/palantir_re_analysis/")
saveRDS(seuratobj,"Cleanupmalatmyoloidwithpseduotime.rds")
FeaturePlot(seuratobj, feature=c("Pseudotime","Entropy"),pt.size = 1)
VlnPlot(seuratobj,features = c("Pseudotime","Entropy"),group.by = "celltype",split.plot = FALSE,split.by = "group",pt.size = 0)
ggsave("Pseudotime_vlnplot_noplot2_splitbygroup2.pdf",width=10,height=5,path = "Fianl analysis_1_19_22/palantir_re_analysis/")
VlnPlot(seuratobj,features = c("Pseudotime","Entropy"),group.by = "celltype",pt.size = 0)
ggsave("Pseudotime_vlnplot_2.pdf",width=5,height=5,path = "Fianl analysis_1_19_22/palantir_re_analysis/")
VlnPlot(seuratobj,features = c("TDC"),group.by = "celltype",split.plot = FALSE,pt.size = 0)
ggsave("DCTerminal state_vlnplot2_4.pdf",width=5,height=5,path = "Fianl analysis_1_19_22/palantir_re_analysis/")

#BOXPLOT Entropy and terminal state probobility
plantir_value<-read.csv("Fianl analysis_1_19_22/palantir_re_analysis/palantir_meta_data_final.csv")
class(plantir_value)
colnames(plantir_value)
library(ggplot2)
ggplot(plantir_value, aes(x = ClusterName, y = Cxcl9.Gbp2b.Mac, fill = Group)) + 
  geom_boxplot() +
  theme(legend.position = "right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
ggsave("Terminstate_box.pdf",width=8,height=5,path = "Fianl analysis_1_19_22/palantir_re_analysis/")



# cluster based pathway analysis
Cluster.markers <- FindAllMarkers(seuratobj, min.pct = 1,only.pos=TRUE)
head(Cluster.markers, n = 20)
write.csv(Cluster.markers,"Cluster.marker_min.pct_1_fianl.csv") 
library(gplots)
library(ggplot2)
sce.markers<-read.csv("Cluster.marker_fianl_fc1.0threshold.csv")
sce.markers$X
library(clusterProfiler)
library(dplyr)
library(ggplot2)
source("kegg_plot_function.R")
#source表示运行整个kegg_plot_function.R脚本，里面是一个function
#以up_kegg和down_kegg为输入数据做图
#symbol transfer to enterzid
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
ids=bitr(sce.markers$gene,'SYMBOL','ENTREZID','org.Mm.eg.db') ## 将SYMBOL转成ENTREZID
sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
View(sce.markers)
saveRDS(sce.markers,"cluster.markerthresholdFC1.0_3_20_22.rds")
## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 
## KEGG
xx <- compareCluster(gcSample,
                     fun = "enrichKEGG",
                     organism = "mmu", pvalueCutoff = 0.05
)
p <- dotplot(xx)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("KEGG_clusterFC1.0.pdf",width = 15,height = 10,path = "Fianl analysis_1_19_22/")
## GO
y <- compareCluster(gcSample,
                    fun = "enrichGO",
                    OrgDb = "org.Mm.eg.db",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.01,
                    qvalueCutoff = 0.05
)
#去除里面含有NA的值否则会引起dotplot报错
y@compareClusterResult<-dplyr::filter(y@compareClusterResult,!is.na(y@compareClusterResult$Description))
#筛选某些Descripition中的某些通路

selected_pathway<-y@compareClusterResult$Description[c(1,2,3,4,192,193,425,427,613,614,665)]
selected_pathway
y@compareClusterResult<-subset(y,y@compareClusterResult$Description %in% selected_pathway)
selected_pathway
write.csv(y@compareClusterResult,"GO_analysis_FC_1.0_des.csv")
p2 <- dotplot(y)
p2 + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))
ggsave("Go_cluster_1.0_selected_final_final.pdf",width = 6,height = 7,path = "Fianl analysis_1_19_22/")
