rm(list = ls())  
load(file = 'step4output.Rdata')
deg<-read.csv("ICI.markers.csv")
deg$symbol<-deg$X
deg$avg_logFC<-deg$avg_log2FC
logFC_t=1
P.Value_t = 0.01
k1 = (deg$p_val_adj < P.Value_t)&(deg$avg_logFC < -logFC_t)
k2 = (deg$p_val_adj < P.Value_t)&(deg$avg_logFC > logFC_t)
change = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change)
table(deg$change)


# volcano plot 

dat<-deg
for_label <- dat%>% 
  filter(symbol %in% c("Cxcl9","Cxcl10","Lgals3","Il1b","Ccrl2","Cbr2","F13a1","Egr1","Pf4","Jun","Lyve1")) 
#%>%
#  head(10)
#for_label <- dat%>% 
#filter(symbol %in% c("Cxcl9","Cxcl10","Lgals3","Il1b","Ccrl2","Cbr2","F13a1","Egr1","Pf4","Jun","Lyve1")) 
#for_label = dat %>% head(10)

x1 = head(dat[dat$change=="up",] ,10)
x2 = tail(dat[dat$change=="down",],4)

for_label = rbind(x1,x2)

dat<-dat[which(dat$p_val_adj != 0),]

P.Value_t=0.01
logFC_t=1
library(ggplot2)

p <- ggplot(data = dat, 
            aes(x = avg_logFC, 
                y = -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
p
volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
volcano_plot+xlim(-3,3)
dat$log<--log10(dat$p_val_adj)
deg2<-dat[order(dat$log),]
ggsave("two group volcano new.png")
write.csv(dat,"DE.ICIvsCtrl.csv")
#将阈值调到0.5 由于差异基因较少
deg<- read.csv("DE.ICIvsCtrl.csv")
logFC_t=0.5
P.Value_t = 0.01
k1 = (deg$p_val_adj < P.Value_t)&(deg$avg_logFC < -logFC_t)
k2 = (deg$p_val_adj < P.Value_t)&(deg$avg_logFC > logFC_t)
change2 = ifelse(k1,"down",ifelse(k2,"up","stable"))
deg <- mutate(deg,change2)
table(deg$change2)

#富集分析考验网速，因此给大家保存了Rdata
#上课运行示例数据无需修改，在做自己的数据时请注意把本行之后的load()去掉
library(clusterProfiler)
library(dplyr)
library(ggplot2)
source("kegg_plot_function.R")
#source表示运行整个kegg_plot_function.R脚本，里面是一个function
#以up_kegg和down_kegg为输入数据做图

#symbol transfer to enterzid

library(org.Mm.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Mm.eg.db)#人类
#其他物种http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))

write.csv(deg,"DE.ICIvsCtrltransfered.csv")




#1.GO database analysis ----

#(1)输入数据
deg<-read.csv("DE.ICIvsCtrltransfered.csv")
gene_up = deg[deg$change2 == 'up','ENTREZID'] 
gene_down = deg[deg$change2 == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#(2)GO分析，分三部分
#以下步骤耗时很长，实际运行时注意把if后面的括号里F改成T
#if(!file.exist("ego_GSE42872.Rdata")) 判断是否存在这个文件
if(T){
  #细胞组分
  ego_CCup <- enrichGO(gene = gene_up, #gene_diff,  or gene_down
                     OrgDb= org.Mm.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  #生物过程
  ego_BPup <- enrichGO(gene = gene_up,#gene_up, #gene_diff,  or gene_down
                     OrgDb= org.Mm.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  #分子功能：
  ego_MFup <- enrichGO(gene = gene_up, #gene_diff,  or gene_down
                     OrgDb= org.Mm.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.01,
                     readable = TRUE)
  save(ego_CCup,ego_BPup,ego_MFup,file = "upego_ICIvsCtrlfcthreshold0.5.Rdata")
}
load(file = "ego_Infvsresfcthreshold0.5.Rdata")

#(3)可视化
#条带图
barplot(ego_MFup,showCategory = 5)
ggsave("barplot_5.png")
#气泡图
dotplot(ego_CCup,x = "GeneRatio")
ggsave("dotplotCCupinres.png",width = 10,height = 8)
dotplot(ego_BPup,x = "GeneRatio")
ggsave("dotplotBPupinres.png",width = 10,height = 8)
dotplot(ego_MFup,x = "GeneRatio")
ggsave("dotplotMFupinres.png",width = 10,height = 8)
library(patchwork)

#下面的图需要映射颜色，设置和示例数据一样的geneList
colnames(deg)
geneList = deg$avg_logFC
names(geneList)=deg$ENTREZID
geneList = sort(geneList,decreasing = T)
#(3)展示top5通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(ego_BPup, showCategory = 10,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
ggsave("ego-MF network.png",width = 15,height = 10)
ggsave("ego-CC network.png",width = 15,height = 10)
ggsave("ego-BPup network.png",width = 15,height = 10)
cnetplot(ego_BPup,showCategory = 10, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
ggsave("ego-BPup network circular.png",width = 20,height = 15)

#Enrichment Map
library(enrichplot)
x2 <- pairwise_termsim(ego_BPup)
emapplot(x2)
#(4)展示通路关系
goplot(ego_BP)
ggsave("pathway network.png")
#(5)Heatmap-like functional classification
heatplot(ego_BPup,foldChange = geneList)
ggsave("ego-BPupheatplot.png",width = 15,height = 10)

#太多基因就会糊。可通过调整比例或者减少基因来控制。
pdf("heatplot.pdf",width = 15,height = 10)
heatplot(ego_CC,foldChange = geneList)
dev.off()

#2.KEGG pathway analysis----
#上调、下调、差异、所有基因
#（1）输入数据
gene_up = deg[deg$change2 == 'up','ENTREZID'] 
gene_down = deg[deg$change2 == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)
gene_all = deg[,'ENTREZID']
#（2）对上调/下调/所有差异基因进行富集分析
#注意这里又有个F
if(T){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'mmu',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'mmu',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'mmu',
                        universe     = gene_all,
                        pvalueCutoff = 0.9)
  save(kk.diff,kk.down,kk.up,file = "GSE42872kegg.Rdata")
}
load("GSE42872kegg.Rdata")
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
g_kegg <- kegg_plot(up_kegg,down_kegg)

g_kegg 
ggsave(g_kegg,filename = 'kegg_up_down.png',width = 20,height = 15)



#gsea作kegg富集分析，可选----
#(1)查看示例数据
data(geneList, package="DOSE")
#(2)将我们的数据转换成示例数据的格式
geneList=deg$avg_logFC
names(geneList)=deg$ENTREZID
geneList=sort(geneList,decreasing = T)
#(3)富集分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'mmu',
                  nPerm        = 1000,
                  minGSSize    = 120,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
#(4)可视化
kegg_plot(up_kegg,down_kegg)
ggsave('kegg_up_down_gsea.png')

