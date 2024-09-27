
load("../TCGA_tumor_tpm.rda")
load("../ssGSEA.rda")

library(WGCNA)
library(tidyverse)

dat <- TCGA[,colnames(TCGA)%in%rownames(mydata)]
PCD_gene <- read.csv("./PCD_1796gene.csv",header = T)
dat <- dat[rownames(dat)%in%PCD_gene$PCD_gene,]

#MAD函数用于剔除表达变化较小的基因，以减少噪音对共表达网络的影响。
datExpr = t(dat[order(apply(dat,1,mad), decreasing = T),])


datExpr <- data.frame(datExpr)
datExpr <- datExpr[complete.cases(datExpr),]
#如果gsg$allOK的结果为TRUE，证明没有缺失值
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK  #TRUE


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK  #TRUE




#聚类所有样本，观察是否有离群值或异常值
sampleTree = hclust(dist(datExpr), method = "average")

par(cex = 0.7);
par(mar = c(0,4,2,0))

pdf(file = '1.Sample clustering to detect outliers.pdf',width = 14,height = 8)

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 0.5,
     cex.axis = 1, cex.main = 1)
cutoff_tree <- 90
# Plot a line to show the cut
abline(h = cutoff_tree, col = "red")
dev.off()

#无离群样本
# Determine cluster under the line  删除离群样本
clust = cutreeStatic(sampleTree, cutHeight = cutoff_tree, minSize = 10)
table(clust) #0：0  1：365
# clust 1 contains the samples we want to keep. 保留clust==1的样本
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)   #1659
nSamples = nrow(datExpr)  #365

#载入表型数据
traitData <- mydata
dim(traitData) #每行是一个样本，每列是一种信息
colnames(traitData) <- "Activity_score"
names(traitData)

#样本信息匹配
traitData <- traitData[rownames(traitData)%in%rownames(datExpr),,drop=F]
table(rownames(datExpr)%in%rownames(traitData))
collectGarbage() #执行垃圾回收，直到空闲内存idicators显示没有更改为止。
#可视化表型数据与基因表达量数据的联系，重构样本聚类树
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(traitData, signed = FALSE) #用颜色代表关联度


pdf(file = '2.Sample dendrogram and trait heatmap.pdf',width = 14,height = 8)

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitData),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
#图片结果解释了临床数据和基因表达量的关联程度，颜色越深，代表这个表型数据与这个样本的基因表达量关系越密切
save(datExpr, traitData, file = "Activity_Liver-01-dataInput.RData")


#构建表达网络
options(stringsAsFactors = FALSE)
enableWGCNAThreads()  #开启多线程
载入第一步保存的数据
lnames = load(file = "Activity_Liver-01-dataInput.RData");
lnames

#构建自动化网络和检测模块
#选择软阈值
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


pdf("3.Threshold.pdf",width = 10, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence")) +
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")+
  abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) +
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

sft$powerEstimate   #合适的软阈值为5 但看图选7
softPower = 7
adjacency = adjacency(datExpr, power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# 调用分层聚类函数
geneTree = hclust(as.dist(dissTOM), method = "average");
# 绘制得出的聚类树（树枝图

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
# 0   1   2   3   4
#245 853 249 186  126 

# Convert numeric lables into colors 将数字标签转换成颜色
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#  grey  turquoise blue brown   yellow
#   245       853   249   186      126   

# Plot the dendrogram and colors underneath

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes 计算模块基因的异质性
MEDiss = 1-cor(MEs);# 计算根据模块特征向量基因计算模块相异度
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf("./4.模块相关系数.pdf")
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      
                      marHeatmap = c(3,4,2,2),
                      
                      plotDendrograms = FALSE,
                      
                      xLabelsAngle = 90)

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


dev.off()


# Call an automatic merging function 调用自动合并功能
# cutHeight 看图筛选想保留的模块 (yellow 变为 brown)
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight =0.3 , verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
pdf(file = "5.geneDendro.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs



#在这个分析中，我们将识别与表型数据显著相关的模块。
#由于我们已经有每个模块的eigengene，只需要将eigengene与外部数据相关联，寻找重要的关联。
datTraits <- traitData
datTraits <- datTraits[rownames(MEs),,drop=F]
# datTraits <- datTraits[,c(1,8,13,26,31)]

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
# 通过相关值对每个关联进行颜色编码
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf("6.Modules.pdf",width = 8,height = 8)
# sizeGrWindow(10,6)
# Will display correlations and their p-values 展示模块与表型数据的相关系数和 P值
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Define variable OS containing the OS column of datTrait
# PI3K = as.data.frame(datTraits$HALLMARK_PI3K_AKT_MTOR_SIGNALING);
# names(PI3K) = "HALLMARK_PI3K_AKT_MTOR_SIGNALING"

Activity = as.data.frame(datTraits$Activity_score);
names(Activity) = 'Activity_score'
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Activity, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Activity), sep="");
names(GSPvalue) = paste("p.GS.", names(Activity), sep="");

# module = "greenyellow"
module = "brown"
# module = "lightyellow"

column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(paste("7.Module Membership in", module,names(Activity), "module.pdf"),width = 8,height = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste0("Gene significance for ",names(Activity)),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#MM-GS图的每一个点：
# 图中的每一个点代表一个基因，横坐标值表示基因与模块的相关性，
# 纵坐标值表示基因与表型性状的相关性，
# 这里可以看出与性状高度显著相关的基因往往是与这个性状显著相关的模块中的重要元素





library(tidyverse)
substanceBXH <- colnames(datExpr) %>% as.character()
gene_symbol <- substanceBXH

annot <- data.frame(substanceBXH = substanceBXH,gene_symbol = gene_symbol )
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       # LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for PI3K
modOrder = order(-abs(cor(MEs, Activity, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")

save(geneInfo,file = 'wgcna_geneInfo.RData')
