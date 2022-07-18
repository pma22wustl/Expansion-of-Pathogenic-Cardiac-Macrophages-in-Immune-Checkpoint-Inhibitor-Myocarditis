rm(list = ls())
library(dplyr)
exp<-read.csv("GSE153981 using data.csv")

# transfer FPKM to TPM


countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#fpkmToTpm
expMatrix<-read.table("ICIM_ICI_DCM_VIM_RPKM.csv",sep = ",",header = T,row.names = 1)
expMatrix2<-expMatrix[!duplicated(expMatrix$SYMBOL),]
expMatrix2<-na.omit(expMatrix2)
rownames(expMatrix2)<-expMatrix2$SYMBOL
expMatrix3<-expMatrix2[,3:31]
expMatrix3<-as.matrix(expMatrix3)
write.csv2(expMatrix3,"FPKM.csv")
class(expMatrix3)
tpms <- apply(expMatrix3,2,fpkmToTpm)

tpms[1:3,]
expMatrix3[1:3,]
colSums(expMatrix3)
#最后可以根据TPM的特征进行检查，看每列加和是否一致
colSums(tpms)
write.csv(tpms_f,"tpm_filtered.csv")
tpms_f=tpms[which(rowSums(tpms) > 0),]
dim(tpms_f)
dim(tpms)
 ex<-tpms
 ex_ICIM<-read.csv("tpm_filtered_ICM.csv",row.names = 1)
#去除行名和列名
#exp2<-exp[,-1]
#rownames(exp2)<-exp2[,1]
#exp3<-exp2[,-1]
#将数据框转换为矩阵
#ex<-as.matrix(exp3)
#将矩阵中字符串转换为数值
#ex=apply(ex,2,as.numeric)
ex=log2(ex+1)
#x<-exp2[,1]
#改matrix ex 的行名,只有矩阵才能做热图
#rownames(ex)<-x
#去除所有值为0的行
ex[ex==0]<-NA
ex=na.omit(ex)##去除所有含NA的行
y<-as.data.frame(ex)
head(ex)
save(ex,y,file="exp for heatmap.Rdata")
exp_final<-y
exp_final$symbol<-row.names(y)
save(exp_final,file="exp_final.Rdata")
library(tidyr)
dat<-gather(data=exp_final,key = sample_nm,value = exp,-symbol)
save(y,dat,file="dat.Rdata")

#subset ICI group
ex_ICI_VM<-ex[,1:13]
ex_ICI_VM<-ex[,-c(14:24)]
ex_ICIM<-ex[]
write.csv(ex_ICI_VM_VIM,"ex_ICI_VM_VIM.csv",sep = ",")
ex_ICI<-read.table("ex_ICI.csv",sep = ",",header = T,row.names = 1)
pl_ICI<-read.csv("Group_list_ICI.csv",sep = ",",header=FALSE,row.names = 1)

# produce grouplist
pl<-read.csv("Group_list.csv",header = TRUE,sep = ",")
colnames(pl_ICI_VM)=c("sample","group")
rownames(pl)<-pl$sample
pl$group
pl_ICI_VM<-pl[-(14:24),]
group_list_ICI_VM<-pl_ICI_VM$group
group_list_ICI<-c(rep("HighCD8",times=5),rep("LowCD8",times=4))
table(pl_ICI_VM$group2)
group_list_ICI
#设置参考水平，对照在前，处理在后
group_list_ICI = factor(group_list_ICI,
                    levels = c("highCD8","LowCD8"))

save(group_list_ICI_VM,file = "group_list_ICI_VM.Rdata")


#设置参考水平，对照在前，处理在后
#转置
group_list_ICI
dat_ICIM=as.data.frame(t(ex_ICIM))

#Boxplot of gene expression
dat_ICIM$group<-group_list_ICI
dat_ICI$group
rownames(dat_ICI)
dat_ICI$CD8A

dat_ICI$CD8<-as.numeric(dat_ICI$CD8)
ggplot(dat_ICI, aes(x = group, y = CD8, fill = group)) + 
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8,
               size = 2, color = "white")
ggsave("2CD8.png",width = 5,height = 3,path = "exported pics/")

#colnames(pl)=c("group","sample")
#rownames(pl)<-pl$sample
#group_list = factor(group_list,
                   # levels = c("Ccl17_CD64low","Ccl17_CD64high"))
#save(exp_final,dat,group_list,pl,file = "Ccl17step1.Rdata")

install.packages("FactoMineR")
library(FactoMineR)#画主成分分析图需要加载这两个包
install.packages("factoextra")
library(factoextra) 
# pca的统一操作走起
dat.pca_ICIM <- PCA(dat_ICIM, graph = FALSE)
pca_plot<-fviz_pca_ind(dat.pca_ICIM,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group_list_ICI, # color by groups
                         palette = c("#00AFBB", "#E7B800"),#修改颜色“”中的一种颜色
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups",repel = "FALSE"
)
pca_plot
ggsave("PCA_ICIMvsICI_2.pdf",width=8,height=5,path = "exported pics/")
save(pca_plot,file = "pca_plot.Rdata")
# PCA on varibles 
dat$group<-group_list
dat[1:3,11891:11892]
dat.pca <- prcomp(dat[,-11892],  scale = TRUE)
fviz_pca_biplot(dat.pca, 
                # Individuals
                geom.ind = "point",
                fill.ind = dat$group, col.ind = "white",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                addEllipses = TRUE,
                # Variables
                alpha.var ="contrib", col.var = "contrib",
                gradient.cols = "RdBu"
)+
  labs(fill = "group", color = "Contrib", alpha = "Contrib") # Change legend title



#热图 
cg=names(tail(sort(apply(ex_ICI,1,sd)),1000))#数据集须用矩阵

n=ex_ICI[cg,]
head(n)

#绘制热图
annotation_col=data.frame(group=group_list_ICI)
rownames(annotation_col)=colnames(n) 
library(pheatmap)
pheatmap(n,
         show_colnames =T,
         show_rownames = F,
         cluster_cols=F,
         annotation_col=annotation_col,
         scale = "row")
library(ggplot2)
ggsave("Heatmap_sd_ICI_highCD8_control.png",width = 15,height = 10,path = "exported pics/")

colnames(pl)=c("group","sample")
rownames(pl)<-pl$sample
#接下来进行3分组的两两之间差异分析
#分析DE genes
library(limma)
design <- model.matrix(~0+group_list_ICI)
colnames(design) <- gsub("group_list_ICI", "", colnames(design))
design
contr.matrix <- makeContrasts(
  ICIMvsICI = ICIM-ICI, 
  levels = colnames(design))
contr.matrix
#生成矩阵
vfit <- lmFit(ex_ICI, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#线性拟合
efit <- eBayes(vfit)
plotSA(efit)
ggsave("efit_ICI_2.png",width=4,height=4,path = "exported pics/")
colnames(efit)
#查看差异结果
summary(decideTests(efit))
save(design,contr.matrix,group_list_ICI,efit,file = "limma different gene expression_ICI_highvscontrol.Rdata")
#挑选两个差异分析结果取交集查看
dt <- decideTests(efit)
summary(dt)
#RT2WvsRT6W和RT2WVSCTRL差异基因的交集
de.common <- which(dt[,1]!=0 & dt[,3]!=0)
length(de.common)
vennDiagram(dt[,c(2,3)], circle.col=c("turquoise", "salmon"))
ggsave("different expressed genes 1.pdf",width=5,height=5,
       path = "exported pics/")
#RT2WVSCTR和RT6WVSCTRL差异基因的交集
de.common_2 <- which(dt[,2]!=0 & dt[,3]!=0)
#统计多少个交集的基因
length(de.common_2)
#Venn图显示
vennDiagram(dt[,2:3], circle.col=c("turquoise", "salmon"))
vennDiagram(dt[,1:3], circle.col=c("turquoise", "salmon"))

#这3次差异分析的结果都是可以独立取出来了：
colnames(efit)
# [1] "AVSB" "AVSC" "BVSC"
ICI_myocarditisvsICI_control <- topTreat(efit, coef=1, n=Inf)
ICI_myocarditisvsVirus_myocarditis <- topTreat(efit, coef=2, n=Inf)
ICI_controlvsVirus_myocarditis <- topTreat(efit, coef=3, n=Inf)
 
save(ICI_myocarditisvsICI_control,
     file = "deg_group_vs_group_ICI2.Rdata")
#coef
#6-9 差异分析
#为deg数据框添加几列 先分析RT6WVSCTRL这一组
deg<-ICI_myocarditisvsICI_control
deg<-ICI_myocarditisvsVirus_myocarditis
deg<-ICI_controlvsVirus_myocarditis
#1.加probe_id列，把行名变成一列
library(dplyr)
deg$symbol <- rownames(deg)

#tibble::rownames_to_column(deg,var="a")
head(deg)
#2.加symbol列，火山图要用
#deg <- inner_join(deg,ids,by="probe_id")
#head(deg)
#按照symbol列去重复
deg <- deg[!duplicated(deg$symbol),]
#3.加change列,标记上下调基因
colnames(deg)
logFC_t=1
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
table(k1)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
change_P.value = ifelse(k1,"down",ifelse(k2,"up","Not Significant"))
deg <- mutate(deg,change_P.value)
table(deg$change_P.value)
deg<-deg[order(deg$logFC),]
colnames(deg)
colnames(deg)[5]="p_val_adj"
write.csv(deg,"ICIM_vsICI.csv")

#1.火山图----
deg<-read.csv("ICI_highHLA-DQA1vsICI_control_3group_ICI.csv",header = T,row.names=1,sep = ",")
dat<-deg
dat<-dat[order(dat$logFC),]
colnames(dat)
colnames(dat)[5]="p_val_adj"
logFC_t=1
P.Value_t = 0.05
n1 = (dat$p_val_adj < P.Value_t)&(dat$logFC < -logFC_t)
n2 = (dat$p_val_adj < P.Value_t)&(dat$logFC > logFC_t)
change = ifelse(n1,"down",ifelse(n2,"up","stable"))
dat <- mutate(dat,change)
table(dat$change)

colnames(dat)
#for_label <- dat%>% 
  #filter(symbol %in% c("MFSD6L","CADPS2","CD19","TSSK6","CCDC13")) 
#%>%
#  head(10)
#for_label <- dat%>% 
#filter(symbol %in% c("HLA-DQA1","Cxcl10","Lgals3","Il1b","Ccrl2","Cbr2","F13a1","Egr1","Pf4","Jun","Lyve1")) 
#for_label = dat %>% head(10)

x1 = tail(dat[dat$change_P.value=="up",] ,99)
write.csv(x1,"ICIMvsICI_UPgenes.csv")
x2 = head(dat[dat$change_P.value=="down",],37)
write.csv(x2,"ICIMvsICIDOWN.csv")
x1 = tail(dat[dat$change=="up",] ,10)
x2 = head(dat[dat$change=="down",],10)
for_label = rbind(x1,x2)
class(for_label)
up_ICIhighvsICI_control<-read.csv("ICI_highCD8vsICI_control_3group_ICI_DOWN.csv")
up_ICIlowvsICI_control<-read.csv("ICI_lowCD8vsICI_control_3group_ICI_DOWN.csv")
up_ICIhighvsICI_low<-read.csv("ICI_highCD8vsICI_loowCD8_3group_ICI_DOWN.csv")
install.packages("VennDiagram")                       # Install VennDiagram package
library("VennDiagram")
#绘制三维韦恩图
venn.plot <- venn.diagram(
  x  = list(
    highvsctrl = up_ICIhighvsICI_control$symbol,
    lowvsctrl = up_ICIlowvsICI_control$symbol,
   highvslow = up_ICIhighvsICI_low$symbol
  ),
  filename = "3triple_Venn_DOWN.png",
  col = "transparent",
  fill = c("red", "blue", "green"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkblue", "white",
                "white", "white", "darkgreen"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 0
)
vennDiagram(up_ICIhighvsICI_control$symbol,up_ICIlowvsICI_control$symbol,
            up_ICIhighvsICI_low$symbol, circle.col=c("turquoise", "salmon"))

UP_common<-intersect(up_ICIhighvsICI_control$X,up_ICIlowvsICI_control$X)
write.csv(UP_common,"high_low_common_genes.csv",sep = ",")
for_label_genes<-read.csv("ICI_highCD8for label.csv",header = FALSE)
for_label_genes[1,1]<-"HLA-DQA1"
for_label<-subset(dat,dat$symbol %in% for_label_genes$V1)
for_label<-dat[dat$symbol=for_label_genes,]
P.Value_t=0.05
logFC_t=1
library(ggplot2)
p <- ggplot(data = dat, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=2.0, 
             aes(color=change_P.value)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
volcano_plot <- p +
  geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    box.padding = 0.3,
    data = for_label,
    color="black",max.overlaps = Inf )
volcano_plot#+xlim(-4,4)+ylim(0,6)

ggsave("Myocarditis vs ealy_phase.pdf",width=8,height=6,path = "exported pics/")

#2.差异基因热图----

#deg<-read.csv("Ccl17_high vs Ccl17_low DE genes.csv")
#colnames(deg)[8]="symbol"
#colnames(deg)[6]="p_val_adj"
#write.csv(deg,"GM vs Sarcoid two group DE genes.csv")
#save(deg,file="Ccl17step2.Rdata")
load("ex deg for heatmap.Rdata")
class(ex)
class(deg)
x=deg$logFC 
names(x)=deg$symbol 
#cg=c(names(head(sort(x),30)),names(tail(sort(x),30)))
#group_list=c(rep("Ccl17_CD64low",times=4),rep("Ccl17_CD64high",times=4))
cg = deg$symbol[deg$change !="stable"]
tail(cg)
genes<-read.csv("log2CPM - for heatmap.csv",header = FALSE)
genes_forheatmap<-genes[3:27,6]
cg<-genes_forheatmap
#select genes in heat map 
cg=c("HLA-DQA1","CXCL10","HLA-DQA1","HLA-DQA1","IFNGR2","HLA-DQA1")
cg1 =c("HLA-DQA1")
for_label_genes_2<-for_label_genes[-(18:22),]
cg1<-for_label_genes_2

cg2<-c(rev(tail(x1$symbol,20)))

#exp_final=exp_final[,order(group_list)]

#group_list_ICI_VM = group_list_ICI_VM[order(group_list_ICI_VM)]
n=ex_ICI[cg2,]
dim(n)
head(n)
save(ex,deg,file="ex deg for heatmap.Rdata")
#作热图
library(pheatmap)
annotation_col=data.frame(group=group_list_ICI)
group_list_ICI
rownames(annotation_col)=colnames(n) 
library(ggplotify)

heatmap_plot <- as.ggplot(pheatmap(n,show_colnames =T,
                                   show_rownames = T,
                                   scale = "row",
                                   cluster_cols = F,
                                   cluster_rows = F,
                                    clustering_method = "average",
                                   annotation_col=annotation_col,
                                   colorRampPalette(c("navy", "white", "firebrick3"))(50)
                                  ,cellwidth = 15, cellheight = 12)) 
# color = colorRampPalette(colors = c("darkblue","white","red"))(100)
heatmap_plot
library("ggplot2")
ggsave(" CD8highvsICIearlyphase_heatmap_20.pdf",width=8,height=6,path = "exported pics/")

load("zero_pca_plot.Rdata")
library(patchwork)
(pca_plot + volcano_plot +heatmap_plot)+ plot_annotation(tag_levels = "A")
ggsave(filename = "Combined.png",width = 25,height = 5,path = "Guoshuai")


#富集分析考验网速，因此给大家保存了Rdata
#上课运行示例数据无需修改，在做自己的数据时请注意把本行之后的load()去掉

library(ggplot2)
library(clusterProfiler)
install.packages("BiocManager")
library(BiocManager)
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
install.packages("clusterProfiler")
library(org.Mm.eg.db)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(bitr)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

write.csv(deg,"DEgens for KeGG.csv")
library(clusterProfiler)
library(dplyr)
library(ggplot2)
source("kegg_plot_function.R")
#source表示运行整个kegg_plot_function.R脚本，里面是一个function
#以up_kegg和down_kegg为输入数据做图

#1.GO database analysis ----

#(1)输入数据
deg<-deg
gene_up = deg[deg$change_P.value == 'up','ENTREZID'] 
gene_down = deg[deg$change_P.value == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#(2)GO分析，分三部分
#以下步骤耗时很长，实际运行时注意把if后面的括号里F改成T
#if(!file.exist("ego_GSE42872.Rdata")) 判断是否存在这个文件
if(T){
  #细胞组分
  ego_CCdown <- enrichGO(gene = gene_down, #gene_diff,  or gene_down
                         OrgDb= org.Hs.eg.db,
                         ont = "CC",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.01,
                         readable = TRUE)
  #生物过程
  ego_BPdown <- enrichGO(gene = gene_down, #gene_diff,  or gene_down
                         OrgDb= org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.01,
                         readable = TRUE)
  #分子功能：
  ego_MFdown <- enrichGO(gene = gene_down, #gene_diff,  or gene_down
                         OrgDb= org.Hs.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         minGSSize = 1,
                         pvalueCutoff = 0.01,
                         qvalueCutoff = 0.01,
                         readable = TRUE)
  save(ego_CCdown,ego_BPdown,ego_MFdown,ego_CCup,ego_BPup,ego_MFup,ego_CCdiff,ego_BPdiff,ego_MFdiff,file = "GO analysis.Rdata")
}
load(file = "GO analysis.Rdata")

#(3)可视化
#条带图
barplot(ego_BPup)
ggsave("barplot_upinSarc.png",path = "Sarc vs GCM")
#气泡图
dotplot(ego_BPup,showCategory=23)
head(ego_BPup@result$Description,30)
result<-ego_BPup@result
write.csv(result,"GO_ego_up.csv",sep = ",")
p1<-dotplot(ego_CCup,x = "GeneRatio")
ggsave("dotplotCCupinCcl17high.png",width = 8,height = 8,path = "GO analysis" )
p2<-dotplot(ego_BPup,x = "GeneRatio")
ggsave("dotplotBPupinCcl17high.png",width = 8,height = 8,path = "GO analysis" )
p3<-dotplot(ego_MFup,x = "GeneRatio")
ggsave("dotplotMFupinSarc.png",width = 8,height = 8,path = "GO analysis" )
library(patchwork)
(p1 + p2 +p3)+ plot_annotation(tag_levels = "A")
ggsave(filename = "Combined.png",width = 30,height = 10,path = "GO analysis")

#下面的图需要映射颜色，设置和示例数据一样的geneList
colnames(deg)
geneList = deg$logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)
#(3)展示top5通路的共同基因，要放大看。
#Gene-Concept Network
ego_BPup
ego_BPup_select<-subset(ego_BPup, ego_BPup@result$Description %in% c("response to interferon-gamma","cellular response to interferon-gamma"))
cnetplot(ego_BPup, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE,circular = FALSE)
ggsave("ego-MF network.png",width = 15,height = 10,path = "GO analysis")
ggsave("ego-CC network.png",width = 15,height = 10,path = "GO analysis")
ggsave("ego-BPup network circle.png",width = 15,height = 10,path = "GO analysis")
cnetplot(ego_BPdown, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
ggsave("ego-BPdown network.png",width = 15,height = 10,path = "GO analysis")

#Enrichment Map
library(enrichplot)
x2 <- pairwise_termsim(ego_BPup)
emapplot(x2)
ggsave("ego-BPup enrichnetwork enrichmap.png",width = 15,height = 10,path = "GO analysis")
#(4)展示通路关系
goplot(ego_BPup)
ggsave("ego-BPdown golplot.png",width = 15,height = 10,path = "GO analysis")
#(5)Heatmap-like functional classification
heatplot(ego_BPup,foldChange = geneList)
ggsave("ego-BPupheatplot.png",width = 15,height = 10)

#太多基因就会糊。可通过调整比例或者减少基因来控制。
pdf("heatplot.pdf",width = 15,height = 10)
heatplot(ego_BPdown,foldChange = geneList)
dev.off()

#2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = deg[deg$change_P.value == 'up','ENTREZID.x'] 
gene_down = deg[deg$change_P.value == 'down','ENTREZID.x'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID.x']
table(deg$change_P.value)
#（2）对上调/下调/所有差异基因进行富集分析
#注意这里又有个F

install.packages("R.utils")

library(R.utils)

R.utils::setOption("clusterProfiler.download.method","auto")


if(T){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9)
  save(kk.diff,kk.down,kk.up,file = "kegg.Rdata")
}
load("kegg.Rdata")
#(3)从富集结果中提取出结果数据框
kegg_diff_dt <- kk.diff@result

#(4)按照pvalue筛选通路
#在enrichkegg时没有设置pvaluecutoff，在此处筛选
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)
#(5)可视化
g_kegg <- kegg_plot(down_kegg,up_kegg)

g_kegg #+ scale_y_continuous(labels = c(15,10,5,0,5,10))
ggsave(g_kegg,filename = 'kegg_up_down.png',width=20,height=20,path = "exported pics/")

up_kegg$Description
up_kegg_select<-up_kegg[c(3,23,35,38,41,45),]
g_kegg_2 <- kegg_plot(up_kegg_select,up_kegg_select)
g_kegg_2
ggsave(g_kegg_2,filename = 'kegg_up_select.png',width=8,height=5,path = "exported pics/")
#gsea作kegg富集分析，可选----
#(1)查看示例数据
data(geneList, package="DOSE")
#(2)将我们的数据转换成示例数据的格式
geneList=deg$logFC
names(geneList)=deg$ENTREZID.x
geneList=sort(geneList,decreasing = T)
#(3)富集分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
save(down_kegg,up_kegg,file="GSEA.Rdata")
#(4)可视化
gse_d<-kegg_plot(up_kegg,down_kegg)
gse_d#+ scale_y_continuous(labels = c(3,2,1,0,1,2,3))
ggsave(filename = 'GSEA.png',width=15,height=10,path = "exported pics/")
