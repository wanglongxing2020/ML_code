library(Seurat)
library(monocle)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggpubr)
library(scales)

#1.读入原始表达数据
HCC1 = Read10X("./dataset/HCC1/")
HCC2 = Read10X("./dataset/HCC2/")
HCC3 = Read10X("./dataset/HCC3/")
# 这里的列名就是barcode
colnames(HCC1) <- paste(colnames(HCC1),"HCC1",sep = "_")
colnames(HCC2) <- paste(colnames(HCC2),"HCC2",sep = "_")
colnames(HCC3) <- paste(colnames(HCC3),"HCC3",sep = "_")

######CCA方法对单细胞数据合并去批次效应######
#数据集中测到的少于200个基因的细胞（min.features = 200）
#和少于3个细胞覆盖的基因（min.cells = 3）被过滤掉
exp1.seurat<- CreateSeuratObject(
  HCC1,
  project = "exp1", 
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
#2.数据质控
exp1.seurat[["percent.mt"]] <- PercentageFeatureSet(exp1.seurat,pattern = "^MT-")
pdf("./plot/2.HCC1_质控1.pdf",width = 8,height = 4)
VlnPlot(exp1.seurat,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
dev.off()

pdf("./plot/2.HCC1_质控2.pdf",width = 8,height = 4)
plot1 <- FeatureScatter(exp1.seurat,feature1 = "nCount_RNA",feature2 = "percent.mt",pt.size = 1.5,group.by = "orig.ident")
plot2 <- FeatureScatter(exp1.seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",pt.size = 1.5,group.by = "orig.ident")
plot1 + plot2
dev.off()



exp1.seurat <- subset(exp1.seurat, 
                      subset = 
                        nFeature_RNA < 3000 & 
                        percent.mt < 10)
dim(exp1.seurat@meta.data)  #15113
#3.标准化
exp1.seurat<-NormalizeData(exp1.seurat,verbose = FALSE)
#鉴定细胞间表达量高变的基因（feature selection）
#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
exp1.seurat<- FindVariableFeatures(exp1.seurat, selection.method = "vst",nfeatures = 2000)

exp2.seurat<- CreateSeuratObject(
  HCC2,
  project = "exp2", 
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
#质控
exp2.seurat[["percent.mt"]] <- PercentageFeatureSet(exp2.seurat,pattern = "^MT-")
pdf("./plot/2.HCC2_质控1.pdf",width = 8,height = 4)
VlnPlot(exp2.seurat,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
dev.off()

pdf("./plot/2.HCC2_质控2.pdf",width = 8,height = 4)
plot1 <- FeatureScatter(exp2.seurat,feature1 = "nCount_RNA",feature2 = "percent.mt",pt.size = 1.5,group.by = "orig.ident")
plot2 <- FeatureScatter(exp2.seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",pt.size = 1.5,group.by = "orig.ident")
plot1 + plot2
dev.off()

exp2.seurat <- subset(exp2.seurat, 
                      subset = 
                        nFeature_RNA < 3000 & 
                        percent.mt < 10)
dim(exp2.seurat@meta.data)  #9060
exp2.seurat<-NormalizeData(exp2.seurat,verbose = FALSE)
exp2.seurat<- FindVariableFeatures(exp2.seurat, selection.method = "vst",nfeatures = 2000)

exp3.seurat<- CreateSeuratObject(
  HCC3,
  project = "exp3", 
  min.cells = 3,
  min.features = 200,
  names.field = 2,
  names.delim = "_")
#质控
exp3.seurat[["percent.mt"]] <- PercentageFeatureSet(exp3.seurat,pattern = "^MT-")
pdf("./plot/2.HCC3_质控1.pdf",width = 8,height = 4)
VlnPlot(exp3.seurat,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
dev.off()

pdf("./plot/2.HCC3_质控2.pdf",width = 8,height = 4)
plot1 <- FeatureScatter(exp3.seurat,feature1 = "nCount_RNA",feature2 = "percent.mt",pt.size = 1.5,group.by = "orig.ident")
plot2 <- FeatureScatter(exp3.seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",pt.size = 1.5,group.by = "orig.ident")
plot1 + plot2
dev.off()

exp3.seurat <- subset(exp3.seurat, 
                      subset = 
                        nFeature_RNA < 3000 & 
                        percent.mt < 10)
dim(exp3.seurat@meta.data)  #14138
exp3.seurat<-NormalizeData(exp3.seurat,verbose = FALSE)
exp3.seurat<- FindVariableFeatures(exp3.seurat, selection.method = "vst",nfeatures = 2000)

#整合
#使用SelectIntegrationFeatures函数选择整合多个样本时使用的特征基因集，通常是每个样本中的高变异基因。
#features <- SelectIntegrationFeatures(object.list = list(exp1.seurat, exp2.seurat, exp3.seurat), nfeatures = 2000)
#author.features = features
exp.anchors<- FindIntegrationAnchors(object.list = c(exp1.seurat,exp2.seurat,exp3.seurat), 
                                     dims = 1:30,
                                     anchor.features=2000)

exp.seurat <- IntegrateData(anchorset = exp.anchors, dims = 1:30)#####包括39311个细胞
save(exp.seurat,file = "exp.seurat.rda")

#4.细胞分类
#（1）分类前首先要对数据集进行降维
#均一化
all.genes <- rownames(exp.seurat)
exp.seurat <- ScaleData(exp.seurat, verbose = FALSE)

pbmc <- exp.seurat
#Perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

pdf("./plot/4.1.1.pdf",width = 12,height = 8)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

pdf("./plot/4.1.2.pdf",width = 12,height = 8)
DimPlot(pbmc, reduction = "pca")
dev.off()

pdf("./plot/4.1.3.pdf",width = 12,height = 8)
DimHeatmap(pbmc, dims = 1, cells = 1000, balanced = TRUE)
dev.off()

pdf("./plot/4.1.4.pdf",width = 12,height = 8)
DimHeatmap(pbmc, dims = 1:15, cells = 1000, balanced = TRUE)
dev.off()


#（2）定义数据集的“维度”
#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pdf("./plot/4.2.PC图.pdf",width = 12,height = 8)
JackStrawPlot(pbmc, dims = 1:15)
dev.off()


#通过筛选贡献小于5%方差和累计贡献90%方差的PC截断点作为曲线拐点
#连续PC之间的差异方差贡献百分比变化小于0.1%的点作为曲线拐点
pct <-pbmc[["pca"]]@stdev /sum(pbmc[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

#选择18
pc.use <-min(which(cumu >90&pct < 5)[1],
             sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)])>0.1),decreasing = T)[1]+1)
pdf("./plot/4.2.碎石图.pdf",width = 12,height = 8)
ElbowPlot(pbmc)$data %>% ggplot()+
  geom_point(aes(x = dims,y=stdev))+
  geom_vline(xintercept = pc.use,color="darkred")+
  theme_bw()+labs(title = "Elbow plot: quantitative approach")
dev.off()


#可以看每个pc的方差
Stdev(pbmc)
save(pbmc,file = "./pbmc.rda")

#（3）细胞分类
pbmc2 <- pbmc
pbmc2 <- FindNeighbors(pbmc, dims = 1:18)
pbmc2 <- FindClusters(pbmc2, resolution = 0.1) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.2) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.3) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.4) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.5) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.6) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.7) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.8) 
pbmc2 <- FindClusters(pbmc2, resolution = 0.9) 
pbmc2 <- FindClusters(pbmc2, resolution = 1) 

##确定分辨率：选择0.3
library(clustree)
library(patchwork)
p1 <- clustree(pbmc2,prefix = "integrated_snn_res.")+coord_flip()
pdf("./plot/4.3.1.pdf",width = 12,height = 10)
p1
dev.off()

#（4）可视化分类结果
#T-SNE
pbmc2 <- RunTSNE(pbmc2, dims = 1:18, label = T)
head(pbmc2@reductions$tsne@cell.embeddings) # 提取T-SNE坐标值。
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
p2 <- DimPlot(pbmc2, reduction = "tsne",group.by = "integrated_snn_res.0.3",label = T,pt.size=0.5)
pdf("./plot/4.4.TSNE.pdf",width = 12,height = 8)
p2
dev.off()

pdf("./plot/4.4.ALL.pdf",width = 18,height = 10)
p1+p2+plot_layout(widths = c(1,1))
dev.off()

saveRDS(pbmc2, file = "pbmc2_tutorial.rds")  #保存数据，用于后续个性化分析

#(5)
#使用ROGUE包,计算单细胞数据中细胞群的纯度（异质性）
#可以看到这里的分群在不同的样本中都具有比较好的均一性
#这里使用“integrated”中的“data”对象
library(ROGUE)
expr <- GetAssayData(pbmc2,layer= "data")%>%as.matrix()
meta <- pbmc2@meta.data
pdf("./plot/4.5.ROGUE.pdf",width = 12,height = 8)
rogue(expr = expr,
      labels = meta$integrated_snn_res.0.3,
      samples = meta$orig.ident,
      platform = "UMI",
      span = 0.6)%>%rogue.boxplot()
dev.off()

#选择分辨率为0.3的分簇
pbmc3 <- readRDS("./pbmc2_tutorial.rds")
pbmc3 <- FindClusters(pbmc3, resolution = 0.3) 


#5.注释
DefaultAssay(pbmc3)<-"RNA"
marker <- read.csv("./LIHC_Marker.csv",header = T)

#差异基因

#注意Layer

pbmc3 <- JoinLayers(pbmc3)

all.markers = FindAllMarkers(pbmc3, only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0.25)

save(all.markers,file = "allmarkers.rda")


cellmarker <-paste0(colnames(marker),"marker")
celltype <-colnames(marker)
for (i in celltype){
  markerlist = toupper(marker[,i])
  markerlist = markerlist[markerlist!='']
  var = paste0(i,'marker')
  assign(var,markerlist)
}

#注释T_cell
DotPlot(pbmc3, features = T_cellmarker)#4,5,6:CD3D CD3E CD3G
#注释Hepatocyte
DotPlot(pbmc3, features = Hepatocytemarker)#13:ALB
#注释B_cell
DotPlot(pbmc3, features = B_cellmarker)#9:MS4A1|CD79A
#注释Endothelial_cell
DotPlot(pbmc3, features = Endothelial_cellmarker)#7:PECAM1
#注释NK_cell
DotPlot(pbmc3, features = NK_cellmarker)#0,1,2,10:NKG7 //GNLY
#注释Macrophage
DotPlot(pbmc3, features = Macrophagemarker)#12:C1QB /C1QC/C1QA
#注释Plasma_cell
DotPlot(pbmc3, features = Plasma_cellmarker)#8:MZB1|IGHG1
#注释Mast_cell
DotPlot(pbmc3, features = Mast_cellmarker)#3:HPGDS
#注释Myeloid dendritic cells    
DotPlot(pbmc3, features = mDCsmarker) #11:CLEC9A

#3、11
cluster <- all.markers[all.markers$cluster=="11",]
cluster <- cluster[order(cluster$avg_log2FC,decreasing = T),]
gene <- cluster$gene[1:10]
gene <- c("NKG7","KLRD1","GNLY")
gene <- "PTGDR2"
DotPlot(pbmc3, features = gene ,
        assay='RNA' ) +
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) +
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #颜色





# 0：NK_cell     1：NK_cell   2：NK_cell   3：Mast_cell    4：T_cell    5： T_cell   
# 6：T_cell   7：Endothelial_cell   8：Plasma_cell 9: B_cell     10: NK_cell
# 11: Myeloid dendritic cell          12: Macrophage          13:Hepatocyte

cluster2celltype = c(
  "0" = "NK cell", 
  "1" = "NK cell", 
  "2" = "NK cell",
  "3" = "Mast cell", 
  "4" = "T cell", 
  "5" = "T cell", 
  "6" = "T cell",
  "7" = "Endothelial cell", 
  "8" = "Plasma cell", 
  "9" = "B cell",
  "10"= "NK cell", 
  "11"= "Myeloid dendritic cell", 
  "12"= "Macrophage", 
  "13"= "Hepatocyte")
pbmc3[['cell_type']] = unname(cluster2celltype[pbmc3@meta.data$seurat_clusters])

pbmc4 <- pbmc3
pbmc4 = RenameIdents(pbmc4,
                           `0`= "NK cell", 
                           `1`= "NK cell", 
                           `2`= "NK cell",
                           `3`= "Mast cell",  
                           `4`= "T cell", 
                           `5`= "T cell", 
                           `6`= "T cell",
                           `7`= "Endothelial cell", 
                           `8`= "Plasma cell", 
                           `9`= "B cell",
                           `10`="NK cell", 
                           `11`="Myeloid dendritic cell", 
                           `12`="Macrophage", 
                           `13`="Hepatocyte")
save(pbmc4,file = 'pbmc4.rda')


DefaultAssay(pbmc4)<-"RNA"
#注释图
pdf('./plot/5.1.细胞注释.pdf',width=12,height = 8)
DimPlot(object = pbmc4, pt.size=0.5,label = T,reduction = "tsne")
dev.off()

#提取默认颜色代码 调整了顺序的
p1 <- DimPlot(object = pbmc4, pt.size=0.5,label = T,reduction = "tsne")
x<-ggplot_build(p1)
info = data.frame(colour = x$data[[1]]$colour, group = x$data[[1]]$group)
info <- unique((arrange(info, group)))
cols <- as.character(info$colour)



#小提琴图

marker <- c("NKG7","HPGDS","CD3D","PECAM1","MZB1","MS4A1","CLEC9A","C1QB","ALB")
pbmc4@meta.data$cell_type <- factor(pbmc4@meta.data$cell_type,levels = rev(c("NK cell","Mast cell","T cell",
                                                                         "Endothelial cell","Plasma cell","B cell",
                                                                         "Myeloid dendritic cell","Macrophage","Hepatocyte")))

pdf('./plot/5.2.标记基因散点图.pdf',width=12,height = 8)

#调整顺序
marker <- c("NKG7","HPGDS","CD3D","PECAM1","MZB1","MS4A1","CLEC9A","C1QB","ALB")
DotPlot(pbmc4, features = marker,group.by = 'cell_type') + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
dev.off()



pdf(file = './plot/5.3.标志基因小提琴图.pdf',width = 12,height = 8) #不常用
VlnPlot(pbmc4, features =marker,pt.size = 0,
        ncol = 3,group.by = 'seurat_clusters')
dev.off()


#标记基因TSNE图
cols2 = c("gray", "coral2")
pdf('./plot/5.4.标记基因TSNE图.pdf',width=12,height = 8)
FeaturePlot(pbmc4,features = marker,cols = cols2,reduction = "tsne")
dev.off()




#占比图
bar.df <- mutate(pbmc4@meta.data)
color_cluster=cols
names(color_cluster)=c("NK cell","Mast cell","T cell",
                       "Endothelial cell","Plasma cell","B cell",
                       "Myeloid dendritic cell","Macrophage","Hepatocyte")
bar.df$orig.ident <- "HCC"
bar.df2 <- mutate(pbmc4@meta.data)

bar.df3 <- rbind(bar.df,bar.df2)



pdf(file = './plot/5.5.proportion.pdf',width = 8,height = 4)
ggplot(bar.df3,aes(x=orig.ident))+
  geom_bar(aes(fill=cell_type),position = "fill",width = .5)+
  scale_x_discrete("")+
  scale_y_continuous("Total cell proportion",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("Cell types",values = color_cluster,limits=names(color_cluster))+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  coord_flip() #让条形图横过来
dev.off()




#三个病人间tSNE图
pdf("./plot/5.6.HCC123.pdf",width = 14,height = 5)
DimPlot(object = pbmc4, pt.size=0.5,label = T,split.by = "orig.ident",reduction = "tsne")
dev.off()


