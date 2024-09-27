library(AUCell)
library(tidyverse)
library(GSEABase)
library(ggraph)
library(Seurat)
#构建gmt文件
rt = read.csv('./mitochondrion.csv',header = T)
gset = c("Hallmark_my_geneset","NA",rt$mitochondrion)
gset = gset%>% 
  as.data.frame() %>% 
  t()
write.table(gset,file = "Hallmark_my_geneset.gmt",sep = "\t",row.names = F,col.names = F,quote = F)

geneSet = getGmt('Hallmark_my_geneset.gmt')
load("~/LIHC_NEW/2.scRNA/pbmc4.rda")

exprMatrix = pbmc4@assays$RNA$counts #RNA 原始值包含全面
exprMatrix = as.matrix(exprMatrix)

#步骤1计算ranking
cells_rankings = AUCell_buildRankings(exprMat = exprMatrix)
#步骤2计算AUC
cell_AUC = AUCell_calcAUC(geneSets = geneSet,cells_rankings)
save(cell_AUC,file = 'cells_AUC.rda')
cell_AUC

#绘图
geneSet = 'Hallmark_my_geneset'
aucs = as.numeric(getAUC(cell_AUC)[geneSet,])
pbmc4$AUC = aucs

df = data.frame(pbmc4@meta.data,pbmc4@reductions$tsne@cell.embeddings)
class_avg = df %>%
  group_by(cell_type) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )

pdf('./AUCell.pdf',width = 12,height = 8)
ggplot(df,aes(tSNE_1,tSNE_2))+
  geom_point(aes(colour = AUC))+
  viridis::scale_color_viridis(option = 'A')+
  ggrepel::geom_label_repel(aes(label = cell_type),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = 'none')+
  theme_bw()
dev.off()

#然后根据中位AUC评分将细胞分为高和低细胞器功能活性 -AUC组
AUC <- df %>%group_by(cell_type)%>%
  summarise(
    AUC_mean = mean(AUC)
  )
median <- median(AUC$AUC_mean)
AUC$group <- ifelse(AUC$AUC_mean < median,"Low","High")
write_csv(AUC,file = "./AUC_mean.csv")

pbmc4$AUC_group <- ifelse(pbmc4$cell_type %in% c("B cell","Mast cell","NK cell","T cell"),"Low","High")
pbmc4$AUC_group <- factor(pbmc4$AUC_group,levels = c("High","Low"))

df1 = data.frame(pbmc4@meta.data,pbmc4@reductions$tsne@cell.embeddings)
class_avg = df1 %>%
  group_by(cell_type) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )


pdf('./AUC_group.pdf',width = 12,height = 8)
ggplot(df1,aes(tSNE_1,tSNE_2))+
  geom_point(aes(colour = AUC_group))+
  scale_color_manual(values = c(High="#FCA311",Low="#E5E5E5"))+
  ggrepel::geom_label_repel(aes(label = cell_type),
                            data = class_avg,
                            size = 6,
                            label.size = 0,
                            segment.color = NA)+
  theme(legend.position = 'none')+
  theme_bw()
dev.off()
pbmc5 <- pbmc4
levels(pbmc5)



 #只能一次RenameIdents
pbmc5 <- readRDS("~/LIHC_NEW/2.scRNA/pbmc2_tutorial.rds")
pbmc5 <- FindClusters(pbmc5, resolution = 0.3) 
DefaultAssay(pbmc5) <- "RNA"
pbmc5 <- JoinLayers(pbmc5)
levels(pbmc5)
pbmc6 = RenameIdents(pbmc5,
                           `0`= "Low", 
                           `1`= "Low", 
                           `2`= "Low",
                           `3`= "Low",  
                           `4`= "Low", 
                           `5`= "Low", 
                           `6`= "Low",
                           `7`= "High", 
                           `8`= "High", 
                           `9`= "Low",
                           `10`="Low", 
                           `11`="High", 
                           `12`="High", 
                           `13`="High")

#differentially active gene  用的slot="data"
DAG.marker = FindAllMarkers(pbmc6, only.pos = TRUE, 
                             min.pct = 0.25, logfc.threshold = 0,assay = "RNA",slot = "data") #1870
DAG.marker2 <- DAG.marker[DAG.marker$p_val_adj < 0.05,] #1711
write.csv(DAG.marker2,file = "./DAG.marker.csv")
