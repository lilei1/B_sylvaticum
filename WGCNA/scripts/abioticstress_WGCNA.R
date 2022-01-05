#!/usr/bin/env Rscript
#   Script to do the WGCNA analysis
# Written by Li Lei, 01-27-2020, Berkley
#getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir);
# Load the WGCNA package
#install.packages("impute")

#install.packages("BiocInstaller")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("impute")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("preprocessCore")

#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("data.table")

#install.packages("WGCNA")
source("https://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)
#dependencies 'robustbase', 'mvtnorm', 'pcaPP' are not available for package 'rrcov'
biocLite("textshape")
biocLite("ggfortify")

install.packages("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel" ) 
install.packages("igraph","ggfortify")
install.packages("ggfortify")
library(igraph)
library(textshape)
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

type = "unsigned"
corType = "pearson"
#corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

tpm <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/normalized_VST_final_shoot_abiotic_formal.txt",sep="\t",header=T)
#head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
metaData <-strsplit(colnames(tpm),"_")
metaData <-data.frame(colnames(tpm))
out <- strsplit(as.character(metaData$colnames.tpm.),'_')
splitted <- data.frame(t(sapply(out, `[`)))
splitted1 <- splitted[,(1-2)]
foo <- data.frame(do.call('rbind', strsplit(as.character(splitted$X3),'.',fixed=TRUE)))
sub_meta <- cbind(splitted$X1,splitted$X2,foo)
colnames(sub_meta) <- c("tissue","treatment","timepoint","replicate")

#metaData <- read.delim("/global/projectb/scratch/llei2019/CSP_Kranthi/Meta_data_susceptible.txt",header = T)
#head(metaData)
#sub_meta <- metaData[(metaData$experiment == "BSMV_MMMV_WSMV"),]
#head(sub_meta)
#nrow(sub_meta)
#col.num <- which(sub_meta$code  %in%  colnames(tpm))
#tpm <- tpm[, col.num]
#ncol(tpm)
#header <- paste(sub_meta$tissue,sub_meta$treatment, sub_meta$timepoint, sub_meta$replicate,sep='_')
#colnames(tpm) <- header
#head(tpm)
femData <- tpm
femData <-data.frame(femData)
### select genes with MAD >85%
dim(femData) #[1] 21794    59
head(femData)
m.mad <- apply(femData,1,mad)#median absolute deviation (MAD)
head(m.mad)
#femData <- femData[which(m.mad > 0),]
#This is to take the 75% genes
#femData <- femData[which(m.mad > 
#above 85%
femData <-  femData[which(m.mad > 
                            max(quantile(m.mad, probs=seq(0, 1, 0.05))[6],0.01)),]
dim(femData)#[1] 16345    59
###
datExpr0 <- as.data.frame(t(femData))
head(datExpr0[,1:8])
dim(datExpr0) #[1]    59 16345
rownames(datExpr0)
colnames(datExpr0)
#Iterations may be required since excluding samples effectively changes criteria on genes and vice versa. 
#The process is repeated until the lists of good samples and genes are stable. 
gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK
### remove the offending genes and samples from the data
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(datExpr0)#[1] 16345
nSamples = nrow(datExpr0)#59
head(datExpr0[,1:8])
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/VST_final_shoot_abiotic_sample_tree_mad75.pdf",width = 25,height = 20)
plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h = 250, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust!=0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)#16345
nSamples = nrow(datExpr)#74

#all of the samples has passed!
#do a PCA analysis for the samples to see what will happen:
ir.pca <- prcomp(datExpr,
                 center = TRUE,
                 scale. = TRUE) 
print(ir.pca)

#head(ir.pca)
head(sub_meta)
sub_meta$group <- paste(sub_meta$tissue,sub_meta$treatment,sub_meta$timepoint, sep='_')

#abiotic_shoot$group <- paste(abiotic_shoot$treatment, abiotic_shoot$timepoint, sep='_')

library(ggfortify)
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/VST_final_shoot_abiotic_pca_mad75.pdf",width = 10,height = 10)
autoplot(ir.pca,data = sub_meta, colour = 'group')
autoplot(ir.pca)
dev.off()
#####
#sampleTree = hclust(dist(datExpr), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
#par(cex = 0.6);
#par(mar = c(0,4,2,0))
#pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_combined_abiotic_distachyon_minor1.pdf",width = 25,height = 20)
#plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#abline(h = 160, col = "red");
#dev.off()


###Automatic block-wise network
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
head(sft)
power = sft$powerEstimate
power
#NA
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/power_soft_thresh_VST_final_shoot_abiotic_mad75.pdf",width = 10,height = 10)
# Scale-free topology fit index as a function of the soft-thresholding power plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.78,col="red")
# this line corresponds to using an R^2 cut-off of h abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#abline(h=0.85,col="red")
dev.off()

#power = 8 is the best according to the plot, I set the mergeCutHeight = 0.3
bwnet = blockwiseModules(datExpr,maxBlockSize = nGenes,
                         power = 7, TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.35, numericLabels = TRUE,pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = "allAin_TOM-blockwise", verbose = 3)
head(bwnet)
bwLabels = bwnet$colors
table(bwnet$colors)
###CutHeight = 0.35 power =8 
#0    1    2    3    4    5    6    7    8    9   10 
#176 6930 2813 2803 2121  793  410  149   75   41   34
#power =7
#0    1    2    3    4    5    6    7    8 
#211 6665 3632 2565 1630 1152  376   74   40 
bwModuleColors = labels2colors(bwLabels)
table(bwModuleColors)
#bwModuleColors #power =8
#black      blue     brown     green      grey   magenta      pink    purple       red 
#149      2813      2803       793       176        41        75        34       410 
#turquoise    yellow 
#6930      2121 
#bwModuleColors #power =7
#black      blue     brown     green      grey      pink       red turquoise    yellow 
#74      3632      2565      1152       211        40       376      6665      1630 
#try to merge:
#merge = mergeCloseModules(datExpr, bwLabels, cutHeight = 0.4, verbose = 3)
#?matchLabels
#moduleLabels = merge$colors
#cbind(moduleLabels, bwnet$colors)
#table(merge$colors)

# open a graphics window
sizeGrWindow(6,6)
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/Dendro_VST_final_shoot_abiotic_mad75_cut0.35_p7.pdf",width = 25,height = 25)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwModuleColors[bwnet$blockGenes[[1]]], mergedColors),
#                    c("Dynamic Tree Cut", "Merged dynamic"), main = "Gene dendrogram and module colors in block 1",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#Power=6 seems to bebetter
#head(bwnet) Since there are not huge different between the merged and unmerged,so still keep the unmerged
#mergedMEs = merge$newMEs
MEs = bwnet$MEs
MEs_col = MEs
#MEs_col = mergedMEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/eigennetwork_VST_final_shoot_abiotic_mad75_cut0.35_p7.pdf",width = 15,height = 13)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()

head(MEs_col)
#graph<-wgcna2igraph(net = bwnet, datExpr = datExpr,
#                    modules2plot = c("blue","green","turquoise","brown"),
#                    colors2plot = c("orange","darkred","cyan","cornflowerblue"),
#                    kME.threshold = 0.5, adjacency.threshold = 0.1,
##                    adj.power = 7, verbose = T,
#                    node.size = 0, frame.color = NA, node.color = NA,
#                    edge.alpha = .5, edge.width =1)
#plot(graph)

#hubs    = chooseTopHubInEachModule(datExpr, 
#                                   bwModuleColors,omitColors = "grey", 
#                                   power = 12, 
#                                 type = "unsigned" )
#This function can only find the single hub genes
#hubs
#trait module association
#In R:
#Plot the relationships between modules and traits
library(dplyr)
#add a stage column for the sub_meta data
header <- paste(sub_meta$tissue,sub_meta$treatment, sub_meta$timepoint, sub_meta$replicate,sep='_')
sub_meta <- mutate(sub_meta, stage = ifelse(sub_meta$timepoint == "1h" | sub_meta$timepoint == "2h", "E", "L"))
sub_meta$group2 <- paste(sub_meta$tissue,sub_meta$treatment,sub_meta$stage,sep='_')
datTraits <- data.frame(header,sub_meta$group2)
str(datTraits)
colnames(datTraits) <- c("ID","group")
str(datTraits)
colnames(tpm)
design = model.matrix(~0+ as.factor(sub_meta$group2))
design <- design[,c(1,3,5,7,2,4,6,8)]
#colnames(design) = c("Bd21_fungi_leaves_Mock","Bd21_fungi_leaves_PCA",
#                     "Bd21-3_bacteria_leaves_Mock","Bd21-3_bacteria_leaves_XT",
#                     "Bd21-3_virus_leaves_Mock","Bd21-3_virus_leaves_BSMV",
#                     "Bd21-3_virus_root_Mock", "Bd21-3_virus_root_BSMV",
#                     "Bd21-3_Endophytes_shoot_Mock","Bd21-3_Endophytes_shoot_SV",
#                     "Bd21-3_Endophytes_root_Mock", "Bd21-3_Endophytes_root_SV")
#colnames(design) = levels(as.factor(datTraits$group))
colnames(design) = c("Shoot_CK_E","Shoot_drought_E","Shoot_heat_E","Shoot_salt_E",
                     "Shoot_CK_L","Shoot_drought_L","Shoot_heat_L", "Shoot_salt_L")
#View(design)
#export the design matrix and see what they look like
#write.table(x = design,
#            file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_beni/design.txt",
#            quote = FALSE,
#            sep = "\t",
#            eol = "\n",
#            col.names = TRUE,
#            row.names = TRUE)

moduleColors = labels2colors(bwnet$colors)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, design, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Display the correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
pdf("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/pheno-module_combined_VST_final_shoot_abiotic_mad75_cut0.35_p7_time_EL.pdf",width = 15,height = 15)
#Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#trait <- "/global/projectb/scratch/llei2019/CSP_Kranthi/pheno_sus_encoding1.txt"

#if(trait != "") {
# traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
#                         check.names=FALSE, comment='',quote="")
# head(traitData)
# sampleName = rownames(datExpr)
# traitData = traitData[match(sampleName, rownames(traitData)), ]
#}
#nrow(MEs_col)
#nrow(traitData)
#robustY = ifelse(corType=="pearson",T,F)
#if (corType=="pearson") {
#  modTraitCor = cor(MEs_col, traitData, use = "p")
#  modTraitP = corPvalueStudent(modTraitCor, nSamples)
#} else {
#  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
#  modTraitCor = modTraitCorP$bicor
#  modTraitP   = modTraitCorP$p
#}

#textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
#dim(textMatrix) = dim(modTraitCor)

#sizeGrWindow(13, 9)
#par(mar = c(5, 30, 5, 5),mfrow = c(1,1))
#cex1 = 1
#pdf("/global/projectb/scratch/llei2019/CSP_Kranthi/pheno-module_combined_abiotic_VST_sus_all_sample_mads75_encoding1_0.25p6.pdf",width = 15,height = 15)
#labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
#               yLabels = colnames(MEs_col), 
#               cex.lab = 1, 
#               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
#               colors = blueWhiteRed(50), 
#               textMatrix = textMatrix, setStdMargins = FALSE, 
#               cex.text = 0.8, zlim = c(-1,1),
#               main = paste("Module-trait relationships"))
#dev.off()

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExpr, MEs, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
head(geneModuleMembership)
head(datExpr[,1:8])
if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(datExpr, design, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExpr, design, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

###Read correlation cyan for heat positive regulation
module = "turquoise"
pheno = "Shoot_drought_L"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/MEturquiose_shoot_drought_L.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
###green
module = "green"
pheno = "Shoot_heat_E"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/MEgreen_shoot_heat_E.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
###brown
module = "brown"
pheno = "Shoot_heat_E_L"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/MEbrown_shoot_heat_E_L.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#####
#####
###Read correlation turquoise for heat positive regulation
#module = "turquoise"
#pheno = "drought"
#modNames = substring(colnames(MEs_col), 3)
#module_column = match(module, modNames)
#pheno_column = match(pheno,colnames(traitData))
#moduleGenes = bwModuleColors == module
#sizeGrWindow(7, 7)
#par(mfrow = c(1,1))
#??geneTraitSignificance
#pdf("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/MEblue_quoise_drought_combined_abiotic_VST_distachyon.pdf",width = 15,height = 15)
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
#                   abs(geneTraitCor[moduleGenes, pheno_column]),
#                   xlab = paste("Module Membership in", module, "module"),
#                   ylab = paste("Gene significance for", pheno),
#                   main = paste("Module membership vs. gene significance\n"),
#                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#dev.off()

###
#In R:
#Intramodular connectivity, module membership, and screening for intramodular hub genes
#Intramodular connectivity
connet=abs(cor(datExpr,use="p"))^7
Alldegrees1=intramodularConnectivity(connet, moduleColors)
head(Alldegrees1)
#Relationship between gene significance and intramodular connectivity
#Here is doing the module and trait association
which.module="turquoise"
drought = as.data.frame(design[,3])
names(drought) = "Shoot_drought"
GS1 = as.numeric(cor(datExpr, drought, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%            1%            2%            3%            4%            5% 
#  -0.9284265518 -0.8425378492 -0.8175020222 -0.7968044441 -0.7804330993 -0.7669607365 
#6%            7%            8%            9%           10%           11% 
#  -0.7539986846 -0.7399321742 -0.7261843035 -0.7150023466 -0.7027678766 -0.6909766290 
#12%           13%           14%           15%           16%           17% 
#  -0.6796069162 -0.6690581879 -0.6573003271 -0.6448939998 -0.6328812518 -0.6225031572 
#18%           19%           20%           21%           22%           23% 
#  -0.6117543590 -0.5991818716 -0.5866124504 -0.5741213826 -0.5603317228 -0.5486811026 
#24%           25%           26%           27%           28%           29% 
#  -0.5366645421 -0.5236716736 -0.5091713115 -0.4959175243 -0.4817607754 -0.4683886670 
#30%           31%           32%           33%           34%           35% 
#  -0.4563510880 -0.4417949751 -0.4258525501 -0.4112529287 -0.3965046359 -0.3820746586 
#36%           37%           38%           39%           40%           41% 
#  -0.3636394818 -0.3472294945 -0.3297682523 -0.3125891096 -0.2950397098 -0.2759176306 
#42%           43%           44%           45%           46%           47% 
#  -0.2577430364 -0.2389101193 -0.2189768169 -0.1996990174 -0.1759070392 -0.1576023765 
#48%           49%           50%           51%           52%           53% 
#  -0.1343024752 -0.1120108243 -0.0921915912 -0.0703389667 -0.0477259660 -0.0234665893 
#54%           55%           56%           57%           58%           59% 
#  -0.0008944749  0.0233015067  0.0485544251  0.0728937862  0.0977536102  0.1226984439 
#60%           61%           62%           63%           64%           65% 
#  0.1495833086  0.1767071206  0.2040239404  0.2285715152  0.2561903637  0.2804720662 
#66%           67%           68%           69%           70%           71% 
#  0.3069776973  0.3304286743  0.3574052940  0.3797754804  0.4059725761  0.4283337194 
#72%           73%           74%           75%           76%           77% 
#  0.4504280004  0.4727083163  0.4928971476  0.5141109688  0.5349107298  0.5529665344 
#78%           79%           80%           81%           82%           83% 
#  0.5725111185  0.5905589003  0.6080541286  0.6225642077  0.6408642179  0.6586029220 
#84%           85%           86%           87%           88%           89% 
#  0.6730812399  0.6890147190  0.7035883600  0.7198681279  0.7354060195  0.7493506277 
#90%           91%           92%           93%           94%           95% 
#  0.7635147981  0.7763932849  0.7904408791  0.8034222314  0.8168366636  0.8324099039 
#96%           97%           98%           99%          100% 
#0.8458546496  0.8603313069  0.8781210475  0.9006741582  0.9578358551

#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
#Display the first few rows of the data frame
head(datKME)
#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
quantile(abs(datKME$MM.turquoise),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#4.181927e-05 1.261057e-02 2.715833e-02 4.234301e-02 5.556366e-02 6.979198e-02 8.246177e-02 
#7%           8%           9%          10%          11%          12%          13% 
#9.587557e-02 1.109221e-01 1.243270e-01 1.360266e-01 1.497717e-01 1.613627e-01 1.730855e-01 
#14%          15%          16%          17%          18%          19%          20% 
#1.857920e-01 1.984012e-01 2.135522e-01 2.253445e-01 2.386473e-01 2.511114e-01 2.632220e-01 
#21%          22%          23%          24%          25%          26%          27% 
#2.760407e-01 2.884476e-01 3.008380e-01 3.126852e-01 3.244429e-01 3.381421e-01 3.496633e-01 
#28%          29%          30%          31%          32%          33%          34% 
#3.613096e-01 3.725305e-01 3.822067e-01 3.952354e-01 4.067334e-01 4.168437e-01 4.272860e-01 
#35%          36%          37%          38%          39%          40%          41% 
#4.372240e-01 4.475030e-01 4.562857e-01 4.671585e-01 4.779974e-01 4.868717e-01 4.971774e-01 
#42%          43%          44%          45%          46%          47%          48% 
#5.062599e-01 5.169094e-01 5.248875e-01 5.354633e-01 5.451847e-01 5.546804e-01 5.636208e-01 
#49%          50%          51%          52%          53%          54%          55% 
#5.726400e-01 5.817488e-01 5.899055e-01 5.986495e-01 6.073780e-01 6.151623e-01 6.234983e-01 
#56%          57%          58%          59%          60%          61%          62% 
#6.308321e-01 6.388921e-01 6.470190e-01 6.542721e-01 6.633243e-01 6.708935e-01 6.791160e-01 
#63%          64%          65%          66%          67%          68%          69% 
#6.865133e-01 6.946350e-01 7.013876e-01 7.094831e-01 7.175326e-01 7.244684e-01 7.304522e-01 
#70%          71%          72%          73%          74%          75%          76% 
#7.394321e-01 7.470693e-01 7.539393e-01 7.612273e-01 7.673274e-01 7.737557e-01 7.799046e-01 
#77%          78%          79%          80%          81%          82%          83% 
#7.853051e-01 7.919381e-01 7.981688e-01 8.048250e-01 8.115517e-01 8.181636e-01 8.257181e-01 
#84%          85%          86%          87%          88%          89%          90% 
#8.318392e-01 8.380331e-01 8.448361e-01 8.510857e-01 8.582375e-01 8.643637e-01 8.706924e-01 
#91%          92%          93%          94%          95%          96%          97% 
#8.779260e-01 8.855334e-01 8.923699e-01 8.996339e-01 9.078128e-01 9.168149e-01 9.274448e-01 
#98%          99%         100% 
#9.379279e-01 9.520147e-01 9.848925e-01 

# if take GS1 99% quantile && datKME$MM.turquoise 99%
FilterGenes= abs(GS1)>0.9 & abs(datKME$MM.turquoise)>0.952
table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#16286    59 
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
#[1] "Brasy1G068600.v1.1" "Brasy1G108000.v1.1" "Brasy1G119300.v1.1" "Brasy1G197800.v1.1"
#[5] "Brasy1G214600.v1.1" "Brasy1G232600.v1.1" "Brasy1G270900.v1.1" "Brasy1G295200.v1.1"
#[9] "Brasy1G357600.v1.1" "Brasy1G456300.v1.1" "Brasy1G529700.v1.1" "Brasy1G530200.v1.1"
#[13] "Brasy1G533800.v1.1" "Brasy1G540900.v1.1" "Brasy2G018800.v1.1" "Brasy2G039700.v1.1"
#[17] "Brasy2G043200.v1.1" "Brasy2G048100.v1.1" "Brasy2G111700.v1.1" "Brasy2G122300.v1.1"
#[21] "Brasy2G129200.v1.1" "Brasy2G137400.v1.1" "Brasy2G187600.v1.1" "Brasy2G247800.v1.1"
#[25] "Brasy2G256300.v1.1" "Brasy2G262600.v1.1" "Brasy2G277100.v1.1" "Brasy2G349700.v1.1"
#[29] "Brasy2G375100.v1.1" "Brasy2G398200.v1.1" "Brasy3G174900.v1.1" "Brasy3G234500.v1.1"
#[33] "Brasy3G247900.v1.1" "Brasy4G043600.v1.1" "Brasy4G110500.v1.1" "Brasy4G199600.v1.1"
#[37] "Brasy4G215400.v1.1" "Brasy4G262500.v1.1" "Brasy4G303500.v1.1" "Brasy5G003000.v1.1"
#[41] "Brasy5G003100.v1.1" "Brasy5G056900.v1.1" "Brasy5G276500.v1.1" "Brasy5G304600.v1.1"
#[45] "Brasy5G398400.v1.1" "Brasy5G514800.v1.1" "Brasy6G081300.v1.1" "Brasy6G268000.v1.1"
#[49] "Brasy7G009200.v1.1" "Brasy8G069700.v1.1" "Brasy8G072700.v1.1" "Brasy8G085600.v1.1"
#[53] "Brasy8G095500.v1.1" "Brasy8G102300.v1.1" "Brasy8G113200.v1.1" "Brasy8G155200.v1.1"
#[57] "Brasy9G084800.v1.1" "Brasy9G125100.v1.1" "Brasy9G193900.v1.1"

#check the control to see if there are any genes overlapped with the control
which.module="turquoise"
CK= as.data.frame(design[,1])
names(CK) = "Shoot_CK"
GS1 = as.numeric(cor(datExpr, CK, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%          1%          2%          3%          4%          5%          6% 
#-0.80208329 -0.61376715 -0.57279975 -0.54859698 -0.52836093 -0.51126814 -0.49504391 
#7%          8%          9%         10%         11%         12%         13% 
#-0.48033262 -0.46753464 -0.45562067 -0.44398101 -0.43530140 -0.42592659 -0.41620036 
#14%         15%         16%         17%         18%         19%         20% 
#-0.40589216 -0.39501543 -0.38587787 -0.37703438 -0.36839228 -0.35844935 -0.35000192 
#21%         22%         23%         24%         25%         26%         27% 
#-0.34129893 -0.33122514 -0.32091503 -0.31069810 -0.29751597 -0.28820319 -0.27699895 
#28%         29%         30%         31%         32%         33%         34% 
#-0.26613386 -0.25470690 -0.24258073 -0.23146630 -0.22001227 -0.20841962 -0.19653720 
#35%         36%         37%         38%         39%         40%         41% 
#-0.18484234 -0.17266180 -0.16012428 -0.14556166 -0.13298783 -0.11912721 -0.10521029 
#42%         43%         44%         45%         46%         47%         48% 
#-0.09271232 -0.07969972 -0.06732038 -0.05382766 -0.04009166 -0.02754492 -0.01442226 
#49%         50%         51%         52%         53%         54%         55% 
#-0.00213612  0.00935005  0.02151486  0.03427922  0.04707446  0.05936843  0.07228535 
#56%         57%         58%         59%         60%         61%         62% 
#0.08232108  0.09332014  0.10445490  0.11643278  0.12716347  0.13796296  0.14717212 
#63%         64%         65%         66%         67%         68%         69% 
#0.15810494  0.16792982  0.18006205  0.19080691  0.20210752  0.21337225  0.22371287 
#70%         71%         72%         73%         74%         75%         76% 
#0.23510167  0.24434847  0.25496793  0.26540861  0.27713896  0.28716298  0.29793384 
#77%         78%         79%         80%         81%         82%         83% 
#0.31029453  0.32140319  0.33112777  0.34243613  0.35181945  0.36170868  0.37342401 
#84%         85%         86%         87%         88%         89%         90% 
#0.38451872  0.39633615  0.40559327  0.41610836  0.42817398  0.44109838  0.45266354 
#91%         92%         93%         94%         95%         96%         97% 
#0.46387957  0.47687991  0.48997353  0.50413632  0.52062748  0.53631979  0.56038875 
#98%         99%        100% 
#0.58806612  0.62979311  0.79302201 
#GS1 98% quantile && datKME$MM.green 98%
FilterGenes= abs(GS1)>0.63 & abs(datKME$MM.turquoise)>0.952
table(FilterGenes)
#FilterGenes
#FALSE 
#16345 
hubgenes <- rownames(datKME)[FilterGenes]
#character(0)
#####conlucion: there is no hub genes which need to excluded from the file

#Brown with heat:
which.module="brown"
heat = as.data.frame(design[,4])
names(heat) = "Shoot_heat"
GS1 = as.numeric(cor(datExpr, heat, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%            1%            2%            3%            4%            5% 
#  -0.9441100677 -0.7989637782 -0.7550643319 -0.7233722423 -0.6934129045 -0.6710068591 
#6%            7%            8%            9%           10%           11% 
#  -0.6506148872 -0.6292283307 -0.6089711066 -0.5883668993 -0.5690397459 -0.5536107089 
#12%           13%           14%           15%           16%           17% 
#  -0.5355590754 -0.5178975414 -0.5030103267 -0.4892479464 -0.4742689804 -0.4605727747 
#18%           19%           20%           21%           22%           23% 
#  -0.4468244924 -0.4321009072 -0.4193386724 -0.4079765429 -0.3948165217 -0.3806659562 
#24%           25%           26%           27%           28%           29% 
#  -0.3664698559 -0.3536632643 -0.3424698484 -0.3301427366 -0.3176663088 -0.3046418129 
#30%           31%           32%           33%           34%           35% 
#  -0.2924205967 -0.2805422648 -0.2678378696 -0.2559042842 -0.2425373754 -0.2306219811 
#36%           37%           38%           39%           40%           41% 
#  -0.2196466689 -0.2075054898 -0.1944429773 -0.1824453598 -0.1702153331 -0.1581023045 
#42%           43%           44%           45%           46%           47% 
#  -0.1454219918 -0.1323252678 -0.1208187033 -0.1070338351 -0.0937744814 -0.0798329625 
#48%           49%           50%           51%           52%           53% 
#  -0.0662529510 -0.0537929658 -0.0406903517 -0.0280700636 -0.0152449615 -0.0007338992 
#54%           55%           56%           57%           58%           59% 
#  0.0140386996  0.0289982095  0.0465937711  0.0617608136  0.0749995177  0.0897085306 
#60%           61%           62%           63%           64%           65% 
#  0.1054378962  0.1199458258  0.1349152789  0.1516398380  0.1686383458  0.1832720833 
#66%           67%           68%           69%           70%           71% 
#  0.1964776764  0.2111556326  0.2263152827  0.2413039294  0.2578598714  0.2750081419 
#72%           73%           74%           75%           76%           77% 
#  0.2906285829  0.3085924486  0.3241957116  0.3408224839  0.3572038386  0.3756797294 
#78%           79%           80%           81%           82%           83% 
#  0.3942703196  0.4088667642  0.4275751914  0.4434409244  0.4598529724  0.4799484066 
#84%           85%           86%           87%           88%           89% 
#  0.4976614328  0.5146543628  0.5325020638  0.5521741496  0.5722147302  0.5946847078 
#90%           91%           92%           93%           94%           95% 
#  0.6150058274  0.6367641026  0.6576480096  0.6787512164  0.7017685890  0.7260906334 
#96%           97%           98%           99%          100% 
#0.7549867254  0.7822808780  0.8178796213  0.8588771075  0.9852195907  
#Generalizing intramodular connectivity for all genes on the array
datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
#Display the first few rows of the data frame
head(datKME)
#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
quantile(abs(datKME$MM.brown),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#2.101729e-05 8.203624e-03 1.960785e-02 2.921679e-02 3.929052e-02 4.806248e-02 5.822280e-02 
#7%           8%           9%          10%          11%          12%          13% 
#6.795341e-02 7.708800e-02 8.544956e-02 9.474840e-02 1.035080e-01 1.138616e-01 1.239887e-01 
#14%          15%          16%          17%          18%          19%          20% 
#1.331480e-01 1.425556e-01 1.506265e-01 1.608801e-01 1.711106e-01 1.805182e-01 1.889854e-01 
#21%          22%          23%          24%          25%          26%          27% 
#1.992504e-01 2.098997e-01 2.184298e-01 2.276060e-01 2.378274e-01 2.468588e-01 2.555072e-01 
#28%          29%          30%          31%          32%          33%          34% 
#2.648667e-01 2.740654e-01 2.831086e-01 2.906224e-01 2.994288e-01 3.085573e-01 3.176051e-01 
#35%          36%          37%          38%          39%          40%          41% 
#3.255625e-01 3.337894e-01 3.425792e-01 3.503228e-01 3.586653e-01 3.675929e-01 3.763499e-01 
#42%          43%          44%          45%          46%          47%          48% 
#3.849989e-01 3.931964e-01 4.009131e-01 4.098476e-01 4.176928e-01 4.262251e-01 4.339115e-01 
#49%          50%          51%          52%          53%          54%          55% 
#4.421111e-01 4.506595e-01 4.580892e-01 4.655287e-01 4.727330e-01 4.802272e-01 4.881828e-01 
#56%          57%          58%          59%          60%          61%          62% 
#4.963749e-01 5.035828e-01 5.113292e-01 5.186790e-01 5.253040e-01 5.333607e-01 5.405489e-01 
#63%          64%          65%          66%          67%          68%          69% 
#5.479009e-01 5.550563e-01 5.625998e-01 5.706165e-01 5.782714e-01 5.861080e-01 5.936798e-01 
#70%          71%          72%          73%          74%          75%          76% 
#6.007201e-01 6.079061e-01 6.151186e-01 6.226772e-01 6.305038e-01 6.374021e-01 6.441312e-01 
#77%          78%          79%          80%          81%          82%          83% 
#6.518017e-01 6.603282e-01 6.675888e-01 6.763741e-01 6.836670e-01 6.927050e-01 7.009585e-01 
#84%          85%          86%          87%          88%          89%          90% 
#7.088462e-01 7.170549e-01 7.266844e-01 7.359142e-01 7.452318e-01 7.552287e-01 7.653140e-01 
#91%          92%          93%          94%          95%          96%          97% 
#7.755807e-01 7.866385e-01 7.972727e-01 8.096760e-01 8.233388e-01 8.377819e-01 8.539984e-01 
#98%          99%         100% 
#8.758075e-01 8.996810e-01 9.773405e-01 

#GS1 99% quantile && datKME$MM.green 99%
FilterGenes= abs(GS1)>0.86 & abs(datKME$MM.brown)>0.9
#GS1 90% quantile && datKME$MM.green 90%
#FilterGenes= abs(GS1)>0.45 & abs(datKME$MM.green)>0.68
table(FilterGenes)
#FilterGenes
#FALSE  TRUE 
#16303    42 
hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
#[1] "Brasy1G077600.v1.1" "Brasy1G251300.v1.1" "Brasy1G258900.v1.1" "Brasy1G527400.v1.1"
#[5] "Brasy1G528000.v1.1" "Brasy2G047300.v1.1" "Brasy2G057200.v1.1" "Brasy2G314600.v1.1"
#[9] "Brasy2G454500.v1.1" "Brasy2G465000.v1.1" "Brasy2G479300.v1.1" "Brasy3G163700.v1.1"
#[13] "Brasy3G295700.v1.1" "Brasy4G034400.v1.1" "Brasy4G073600.v1.1" "Brasy4G166700.v1.1"
#[17] "Brasy4G168000.v1.1" "Brasy4G201100.v1.1" "Brasy4G389900.v1.1" "Brasy5G055400.v1.1"
#[21] "Brasy5G079000.v1.1" "Brasy5G085200.v1.1" "Brasy5G145300.v1.1" "Brasy5G161200.v1.1"
#[25] "Brasy5G462500.v1.1" "Brasy6G044800.v1.1" "Brasy6G182200.v1.1" "Brasy6G233700.v1.1"
#[29] "Brasy7G178800.v1.1" "Brasy8G124300.v1.1" "Brasy8G174800.v1.1" "Brasy8G228500.v1.1"
#[33] "Brasy9G106700.v1.1" "Brasy9G130100.v1.1" "Brasy9G171300.v1.1" "Brasy9G270400.v1.1"
#[37] "Brasy9G313200.v1.1" "Brasy9G336400.v1.1" "Brasy9G360300.v1.1" "Brasy9G362900.v1.1"
#[41] "BrasyJ025700.v1.1"  "BrasyJ089800.v1.1" 

#check the control
which.module="brown"
heatCK= as.data.frame(design[,1])
names(heatCK) = "Shoot_CK"
GS1 = as.numeric(cor(datExpr, heatCK, use = "p"))
GeneSignificance=abs(GS1)
quantile(GS1,probs = seq(0, 1, 0.01))
#0%          1%          2%          3%          4%          5%          6% 
#-0.80208329 -0.61376715 -0.57279975 -0.54859698 -0.52836093 -0.51126814 -0.49504391 
#7%          8%          9%         10%         11%         12%         13% 
#-0.48033262 -0.46753464 -0.45562067 -0.44398101 -0.43530140 -0.42592659 -0.41620036 
#14%         15%         16%         17%         18%         19%         20% 
#-0.40589216 -0.39501543 -0.38587787 -0.37703438 -0.36839228 -0.35844935 -0.35000192 
#21%         22%         23%         24%         25%         26%         27% 
#-0.34129893 -0.33122514 -0.32091503 -0.31069810 -0.29751597 -0.28820319 -0.27699895 
#28%         29%         30%         31%         32%         33%         34% 
#-0.26613386 -0.25470690 -0.24258073 -0.23146630 -0.22001227 -0.20841962 -0.19653720 
#35%         36%         37%         38%         39%         40%         41% 
#-0.18484234 -0.17266180 -0.16012428 -0.14556166 -0.13298783 -0.11912721 -0.10521029 
#42%         43%         44%         45%         46%         47%         48% 
#-0.09271232 -0.07969972 -0.06732038 -0.05382766 -0.04009166 -0.02754492 -0.01442226 
#49%         50%         51%         52%         53%         54%         55% 
#-0.00213612  0.00935005  0.02151486  0.03427922  0.04707446  0.05936843  0.07228535 
#56%         57%         58%         59%         60%         61%         62% 
#0.08232108  0.09332014  0.10445490  0.11643278  0.12716347  0.13796296  0.14717212 
#63%         64%         65%         66%         67%         68%         69% 
#0.15810494  0.16792982  0.18006205  0.19080691  0.20210752  0.21337225  0.22371287 
#70%         71%         72%         73%         74%         75%         76% 
#0.23510167  0.24434847  0.25496793  0.26540861  0.27713896  0.28716298  0.29793384 
#77%         78%         79%         80%         81%         82%         83% 
#0.31029453  0.32140319  0.33112777  0.34243613  0.35181945  0.36170868  0.37342401 
#84%         85%         86%         87%         88%         89%         90% 
#0.38451872  0.39633615  0.40559327  0.41610836  0.42817398  0.44109838  0.45266354 
#91%         92%         93%         94%         95%         96%         97% 
#0.46387957  0.47687991  0.48997353  0.50413632  0.52062748  0.53631979  0.56038875 
#98%         99%        100% 
#0.58806612  0.62979311  0.79302201 

#GS1 99% quantile && datKME$MM.green 98%
FilterGenes= abs(GS1)>0.63 & abs(datKME$MM.brown)>0.9
table(FilterGenes)
#FilterGenes
#FALSE 
#16345 
####so all of the genes listed are the hub genes

#In R:
#Export the network
#Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 7) 
# Select module
module = "turquoise" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])

#Screen the top genes
nTop = 59
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]

cyt = exportNetworkToCytoscape(filter,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-filter-top59-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-filter-top59-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = rownames(filter), 
                               nodeAttr = moduleColors[inModule][1:nTop])

# Brown
module = "brown" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])

#Screen the top genes
nTop = 42
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]

cyt = exportNetworkToCytoscape(filter,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-filter-top42-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-filter-top42-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = rownames(filter), 
                               nodeAttr = moduleColors[inModule][1:nTop])
#####yellow
module = "yellow" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])
###
#####blue
module = "blue" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])
#####pink
module = "pink" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])
#####red
module = "red" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])
#####black
module = "black" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])
#####green
module = "green" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])

#####grey
module = "grey" #change to green
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule] 
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

#Export the network into edge and node list files Cytoscape can read
#default threshold = 0.5, we could adjust parameter based on actual situations or in Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes, 
                               nodeAttr = moduleColors[inModule])
#####extract top 50 hub genes for each module!!!
MEs <- moduleEigengenes(datExpr, bwModuleColors)$eigengenes
head(MEs)
kMEs <- signedKME(datExpr, MEs)
head(kMEs)
# rank the genes for each module on kMEs
rankGenes <- function(x){
  kMErank <- rank(-kMEs[ ,x])
  genes <- rownames(kMEs)
  genes <- genes[order(kMErank)]
  #genes[1:50]#top 50 hub genes!!!
}

topGenes <- lapply(1:ncol(kMEs), rankGenes)
nrow(topGenes)
# Get the top results in a data.frame
topGenes <- do.call(cbind, topGenes)
colnames(topGenes) <- substr(colnames(kMEs), start=4, stop=30)
nrow(topGenes)
head(topGenes)
hubgenes <- data.frame(topGenes)
write.csv(hubgenes, file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/gene_membership_shoot_abioticstress_mads75_0.35power7.csv")

#
#adj_mat<-adjacency(datExpr,power=12)
adj <- TOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
V(network)$color <- bwnet$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
plot(network, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

###Plot the network for all of the modules
top.n.edges = 25000#3000#5000#15000#25000
min.edge = 2
adj_mat<-adjacency(datExpr,power=7)
head(adj_mat)
nrow(adj_mat)
message("number of genes present = ", nrow(adj_mat))#number of genes present = 16345 
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)
#adjacency.threshold = 0.893484769264947#if edge=5000,0.879482033215977；edge=10000，adjacency.threshold = 0.856292297792383,edge15000,adjacency.threshold = 0.839325416948576
#25000，adjacency.threshold = 0.67461091429808,20000,adjacency.threshold = 0.826031854913583
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#1975, 25000; #356 #461#661#839#1064#962
test <- data.frame(adj_mat)
nrow(test)#356#461#661#839#1064
#row.names(test)[modProbes,]

#turquoise
probes = colnames(datExpr)
turquoise_modules = c("turquoise")
turquoise_inModule = is.finite(match(bwModuleColors, turquoise_modules))
turquoise_modProbes = probes[turquoise_inModule]
turquoise_list <- subset(adj_mat, rownames(adj_mat) %in% turquoise_modProbes)

#pink
probes = colnames(datExpr)
pink_modules = c("pink")
pink_inModule = is.finite(match(bwModuleColors, pink_modules))
pink_modProbes = probes[pink_inModule]
pink_list <- subset(adj_mat, rownames(adj_mat) %in% pink_modProbes)#0

#brown
probes = colnames(datExpr)
brown_modules = c("brown")
brown_inModule = is.finite(match(bwModuleColors, brown_modules))
brown_modProbes = probes[brown_inModule]
brown_list <- subset(adj_mat, rownames(adj_mat) %in% brown_modProbes)

#red
probes = colnames(datExpr)
red_modules = c("red")
red_inModule = is.finite(match(bwModuleColors, red_modules))
red_modProbes = probes[red_inModule]
red_list <- subset(adj_mat, rownames(adj_mat) %in% red_modProbes)

#magenta
#probes = colnames(datExpr)
#magenta_modules = c("magenta")
#magenta_inModule = is.finite(match(bwModuleColors, magenta_modules))
#magenta_modProbes = probes[magenta_inModule]
#magenta_list <- subset(adj_mat, rownames(adj_mat) %in% magenta_modProbes)#0

#blue
probes = colnames(datExpr)
blue_modules = c("blue")
blue_inModule = is.finite(match(bwModuleColors, blue_modules))
blue_modProbes = probes[blue_inModule]
blue_list <- subset(adj_mat, rownames(adj_mat) %in% blue_modProbes)
#green
probes = colnames(datExpr)
green_modules = c("green")
green_inModule = is.finite(match(bwModuleColors, green_modules))
green_modProbes = probes[green_inModule]
green_list <- subset(adj_mat, rownames(adj_mat) %in% green_modProbes)

#brown
#probes = colnames(datExpr)
#brown_modules = c("brown")
#brown_inModule = is.finite(match(bwModuleColors, brown_modules))
#brown_modProbes = probes[brown_inModule]
#brown_list <- subset(adj_mat, rownames(adj_mat) %in% brown_modProbes)

#yellow
probes = colnames(datExpr)
yellow_modules = c("yellow")
yellow_inModule = is.finite(match(bwModuleColors, yellow_modules))
yellow_modProbes = probes[yellow_inModule]
yellow_list <- subset(adj_mat, rownames(adj_mat) %in% yellow_modProbes)

#black
probes = colnames(datExpr)
black_modules = c("black")
black_inModule = is.finite(match(bwModuleColors, black_modules))
black_modProbes = probes[black_inModule]
black_list <- subset(adj_mat, rownames(adj_mat) %in% black_modProbes)#0

#greenyellow
#probes = colnames(datExpr)
#greeny_modules = c("greenyellow")
#greeny_inModule = is.finite(match(bwModuleColors, greeny_modules))
#greeny_modProbes = probes[greeny_inModule]
#greeny_list <- subset(adj_mat, rownames(adj_mat) %in% greeny_modProbes)#0

#grey
probes = colnames(datExpr)
grey_modules = c("grey")
grey_inModule = is.finite(match(bwModuleColors, grey_modules))
grey_modProbes = probes[grey_inModule]
grey_list <- subset(adj_mat, rownames(adj_mat) %in% grey_modProbes)#0

#lightcyan
#probes = colnames(datExpr)
#lightcyan_modules = c("lightcyan")
#lightcyan_inModule = is.finite(match(bwModuleColors, lightcyan_modules))
#lightcyan_modProbes = probes[lightcyan_inModule]
#lightcyan_list <- subset(adj_mat, rownames(adj_mat) %in% lightcyan_modProbes)#0


#midnightblue
#probes = colnames(datExpr)
#midnightblue_modules = c("midnightblue")
#midnightblue_inModule = is.finite(match(bwModuleColors, midnightblue_modules))
#midnightblue_modProbes = probes[midnightblue_inModule]
#midnightblue_list <- subset(adj_mat, rownames(adj_mat) %in% midnightblue_modProbes)#0



#salmon
#probes = colnames(datExpr)
#salmon_modules = c("salmon")
#salmon_inModule = is.finite(match(bwModuleColors, salmon_modules))
#salmon_modProbes = probes[salmon_inModule]
#salmon_list <- subset(adj_mat, rownames(adj_mat) %in% salmon_modProbes)#0

#tan
#probes = colnames(datExpr)
#tan_modules = c("tan")
#tan_inModule = is.finite(match(bwModuleColors, tan_modules))
#tan_modProbes = probes[tan_inModule]
#tan_list <- subset(adj_mat, rownames(adj_mat) %in% tan_modProbes)#0



#darkgreen
#probes = colnames(datExpr)
#darkgreen_modules = c("darkgreen")
#darkgreen_inModule = is.finite(match(bwModuleColors, darkgreen_modules))
#darkgreen_modProbes = probes[darkgreen_inModule]
#darkgreen_list <- subset(adj_mat, rownames(adj_mat) %in% darkgreen_modProbes)
#darkgrey
#probes = colnames(datExpr)
#darkgrey_modules = c("darkgrey")
#darkgrey_inModule = is.finite(match(bwModuleColors, darkgrey_modules))
#darkgrey_modProbes = probes[darkgrey_inModule]
#darkgrey_list <- subset(adj_mat, rownames(adj_mat) %in% darkgrey_modProbes)
#darkred
#probes = colnames(datExpr)
#darkred_modules = c("darkred")
#darkred_inModule = is.finite(match(bwModuleColors, darkred_modules))
#darkred_modProbes = probes[darkred_inModule]
#darkred_list <- subset(adj_mat, rownames(adj_mat) %in% darkred_modProbes)
#darkturquoise
#probes = colnames(datExpr)
#darkturquoise_modules = c("darkturquoise")
#darkturquoise_inModule = is.finite(match(bwModuleColors, darkturquoise_modules))
#darkturquoise_modProbes = probes[darkturquoise_inModule]
#darkturquoise_list <- subset(adj_mat, rownames(adj_mat) %in% darkturquoise_modProbes)
#grey60
#probes = colnames(datExpr)
#grey60_modules = c("grey60")
#grey60_inModule = is.finite(match(bwModuleColors, grey60_modules))
#grey60_modProbes = probes[grey60_inModule]
#grey60_list <- subset(adj_mat, rownames(adj_mat) %in% grey60_modProbes)
#lightgreen
#probes = colnames(datExpr)
#lightgreen_modules = c("lightgreen")
#lightgreen_inModule = is.finite(match(bwModuleColors, lightgreen_modules))
#lightgreen_modProbes = probes[lightgreen_inModule]
#lightgreen_list <- subset(adj_mat, rownames(adj_mat) %in% lightgreen_modProbes)
#lightyellow
#probes = colnames(datExpr)
#lightyellow_modules = c("lightyellow")
#lightyellow_inModule = is.finite(match(bwModuleColors, lightyellow_modules))
#lightyellow_modProbes = probes[lightyellow_inModule]
#lightyellow_list <- subset(adj_mat, rownames(adj_mat) %in% lightyellow_modProbes)

#orange
#probes = colnames(datExpr)
#orange_modules = c("orange")
#orange_inModule = is.finite(match(bwModuleColors, orange_modules))
#orange_modProbes = probes[orange_inModule]
#orange_list <- subset(adj_mat, rownames(adj_mat) %in% orange_modProbes)
#royalblue
#probes = colnames(datExpr)
#royalblue_modules = c("royalblue")
#royalblue_inModule = is.finite(match(bwModuleColors, royalblue_modules))
#royalblue_modProbes = probes[royalblue_inModule]
#royalblue_list <- subset(adj_mat, rownames(adj_mat) %in% royalblue_modProbes)

#rownames(yellow_list) <- subset(yellow_list, colnames(yellow_list) %in% modProbes)
#nrow(yellow_list)
#yellow_list$color <- rep("yellow", nrow(yellow_list))
#yellow_sub <- data.frame(row.names(yellow_list),yellow_list$color)
#network <- graph.adjacency(yellow_list)
#head(network)
#network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
#par(mar=c(0,0,0,0))
# remove unconnected nodes
#network <- delete.vertices(network, degree(network)==0)
#pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network, vertex.color="yellow",vertex.label="",vertex.size=2, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

#dev.off()
####
network <- graph.adjacency(adj_mat,mode = "undirected",weighted = TRUE,)
head(network)
E(network)
V(network)[rownames(turquoise_list)]$color <- "turquoise"
V(network)[rownames(pink_list)]$color <- "pink"
#V(network)[rownames(purple_list)]$color <- "purple"
V(network)[rownames(red_list)]$color <- "red"
#V(network)[rownames(magenta_list)]$color <- "magenta"
V(network)[rownames(blue_list)]$color <- "blue"
V(network)[rownames(green_list)]$color <- "green"
V(network)[rownames(brown_list)]$color <- "brown"
V(network)[rownames(yellow_list)]$color <- "yellow"
V(network)[rownames(black_list)]$color <- "black"
#V(network)[rownames(greeny_list)]$color <- "greenyellow"
V(network)[rownames(grey)]$color <- "grey"
#V(network)[rownames(cyan_list)]$color <- "cyan"
#V(network)[rownames(grey60_list)]$color <- "grey60"
#V(network)[rownames(tan_list)]$color <- "tan"
#V(network)[rownames(midnightblue_list)]$color <- "midnightblue"
#V(network)[rownames(lightcyan_list)]$color <- "lightcyan"
#V(network)[rownames(salmon_list)]$color <- "salmon"
#V(network)[rownames(darkgreen_list)]$color <- "darkgreen"
#V(network)[rownames(darkgrey_list)]$color <- "darkgrey"
#V(network)[rownames(darkred_list)]$color <- "darkred"
#V(network)[rownames(darkturquoise_list)]$color <- "darkturquoise"
#V(network)[rownames(lightgreen_list)]$color <- "lightgreen"
#V(network)[rownames(lightyellow_list)]$color <- "lightyellow"
#V(network)[rownames(orange_list)]$color <- "orange"
#V(network)[rownames(royalblue_list)]$color <- "royalblue"
#all[rownames(yellow_list)]
#all[rownames(yellow_list)]
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
sizeGrWindow(10,10)
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
#minC <- rep(-Inf, vcount(network))
minC <- rep(-Inf, vcount(network))
maxC <- rep(Inf, vcount(network))
minC[1] <- maxC[1] <- 0
#co <- layout_with_fr(network, minx=minC, maxx=maxC,
#                     miny=minC, maxy=maxC)
co <- layout.fruchterman.reingold(network, minx=minC, maxx=maxC,
                                  miny=minC, maxy=maxC)
#layout.fruchterman.reingold
co[1,]
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/network_syl_shoot_abiotic_dges25000.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network,layout=co,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network,layout=layout.graphopt,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5,edge.arrow.size = 0.2)

dev.off()

#Since the turquoise (Drought) and cyan (heat) associated with the certain abiotic treatment,
#so I have to do the GO analysis and the plot the network for those two modules

###Find all of the genes in the module:
gs<-colnames(datExpr)
#  cols<-net[[1]]
cols <- bwModuleColors
names(cols)<-gs #assign each gene into different module
head(gs)

kme <- signedKME(datExpr = datExpr, datME = MEs, outputColumnName = NULL)
head(kme)

#row.names(kme)<-gs
kmes<-sapply(1:length(gs), function(x) abs(kme[gs[x],colnames(kme) == cols[[x]]]))
kmes<-data.frame(genes = gs, cols = cols, kme = kmes)
head(kmes)
nrow(kmes)
#check the gene list for each module
#turquoise
turquoise_kmes <- kmes[kmes$cols=="turquoise",]
head(turquoise_kmes)
nrow(turquoise_kmes)
#6244
turquoise_kmes <-turquoise_kmes[order(-turquoise_kmes$kme),]
write.csv(turquoise_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_turquoise_kmes.csv")
###pink
pink_kmes <- kmes[kmes$cols=="pink",]
head(pink_kmes)
nrow(pink_kmes)
#347
pink_kmes <-pink_kmes[order(-pink_kmes$kme),]
write.csv(pink_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_pink_kmes.csv")

###purple
purple_kmes <- kmes[kmes$cols=="purple",]
head(purple_kmes)
nrow(purple_kmes)
#178
purple_kmes <-purple_kmes[order(-purple_kmes$kme),]
write.csv(purple_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_purple_kmes.csv")

###red
red_kmes <- kmes[kmes$cols=="red",]
head(red_kmes)
nrow(red_kmes)
#444
red_kmes <-red_kmes[order(-red_kmes$kme),]
write.csv(red_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_red_kmes.csv")

###Magenta
magenta_kmes <- kmes[kmes$cols=="magenta",]
head(magenta_kmes)
nrow(magenta_kmes)
#265
magenta_kmes <-magenta_kmes[order(-magenta_kmes$kme),]
write.csv(magenta_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_magenta_kmes.csv")

###Blue
blue_kmes <- kmes[kmes$cols=="blue",]
head(blue_kmes)
nrow(blue_kmes)
#4203
blue_kmes <-blue_kmes[order(-blue_kmes$kme),]
write.csv(blue_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_blue_kmes.csv")

###green
green_kmes <- kmes[kmes$cols=="green",]
head(green_kmes)
nrow(green_kmes)
#587
green_kmes <-green_kmes[order(-green_kmes$kme),]
write.csv(green_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_green_kmes.csv")

###brown
brown_kmes <- kmes[kmes$cols=="brown",]
head(brown_kmes)
nrow(brown_kmes)
#3801
brown_kmes <-brown_kmes[order(-brown_kmes$kme),]
write.csv(brown_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_brown_kmes.csv")

###yellow
yellow_kmes <- kmes[kmes$cols=="yellow",]
head(yellow_kmes)
nrow(yellow_kmes)
#2316
yellow_kmes <-yellow_kmes[order(-yellow_kmes$kme),]
write.csv(yellow_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_yellow_kmes.csv")

###black
black_kmes <- kmes[kmes$cols=="black",]
head(black_kmes)
nrow(black_kmes)
#354
black_kmes <-black_kmes[order(-black_kmes$kme),]
write.csv(black_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_black_kmes.csv")

###greenyellow
greenyellow_kmes <- kmes[kmes$cols=="greenyellow",]
head(greenyellow_kmes)
nrow(greenyellow_kmes)
#52
greenyellow_kmes <-greenyellow_kmes[order(-greenyellow_kmes$kme),]
write.csv(greenyellow_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_greenyellow_kmes.csv")

###grey
grey_kmes <- kmes[kmes$cols=="grey",]
head(grey_kmes)
nrow(grey_kmes)
#100
grey_kmes <-grey_kmes[order(-grey_kmes$kme),]
write.csv(grey_kmes, file = "/global/projectb/scratch/llei2019/CSP_Kranthi/sus_brachy_grey_kmes.csv")

#################
turquoise_gene_list <- rownames(turquoise_kmes)

turquoise_kmes <-turquoise_kmes[order(-turquoise_kmes$kme),]
turquoise_top20_gene_list <- turquoise_kmes[1:20,]$genes

gs<-kmes$genes[kmes$cols=="turquoise"]
#cols<-cols[kmes$kme>=kME.threshold]
datExpr_turquoise = datExpr[,gs]
head(datExpr_turquoise[,1:8])
ncol(datExpr_turquoise)
turquoise_gene_list <- data.frame(colnames(datExpr_turquoise))
#}#yellow module we have 3587 genes!!! 3951 genes for turquoise
###write into the file:
write.csv(turquoise_gene_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_turtoise_list_VST_distachyon.csv")

adj_mat<-adjacency(datExpr_turquoise,power=9)
###
head(adj_mat[,1:8])
######### Get the top edges of the blue module:
top.n.edges = 3000
min.edge = 2
#head(adj_mat)
message("number of genes present = ", nrow(adj_mat))#number of genes present = 3587 Distachyon: 3951
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)#adjacency.threshold = 0.690966673259787#distachyon:0.511680160681778
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#number of genes present = 343
#distachyon: 351
#test <- data.frame(adj_mat)
#row.names(test)[modProbes,])
#yellow
#yellow_tophub <- subset(adj_mat, rownames(adj_mat) %in% yellow_top20_gene_list)
#`%notin%` <- Negate(`%in%`)#build a function not in
#yellow_NOTtophub <- subset(adj_mat, rownames(adj_mat) %notin% yellow_top20_gene_list)

#row.names(yellow_tophub)
#nrow(yellow_tophub)
#yellow_tophub_size <- cbind(row.names(yellow_tophub),rep(7,nrow(yellow_tophub)))
#yellow_NOTtophub_size <- cbind(row.names(yellow_NOTtophub),rep(2.5,nrow(yellow_NOTtophub)))
#yellow_size <-rbind(yellow_tophub_size,yellow_NOTtophub_size)
#yellow_tophub_size<-setNames(row.names(yellow_tophub),rep(7,nrow(yellow_tophub)))

#adj_mat
network <- graph.adjacency(adj_mat,mode = "undirected",weighted = TRUE)
head(network)
E(network)
V(network)
deg <- degree(network, mode="all")
V(network)$vertex_degree <-  deg

network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
V(network)$color <- "turquoise"
scale_factor <- 0.1
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#_turtoise.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=V(network)$vertex_degree*scale_factor, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
dev.off()

#cyan
gs<-colnames(datExpr)
cols <- bwModuleColors
names(cols)<-gs #assign each gene into different module

kme <- signedKME(datExpr = datExpr, datME = MEs, outputColumnName = NULL)
head(kme)
kmes<-sapply(1:length(gs), function(x) abs(kme[gs[x],colnames(kme) == cols[[x]]]))
kmes<-data.frame(genes = gs, cols = cols, kme = kmes)
head(kmes)
cyan_kmes <- kmes[kmes$cols=="cyan",]
head(cyan_kmes)
cyan_kmes <-cyan_kmes[order(-cyan_kmes$kme),]
cyan_top20_gene_list <- cyan_kmes[1:20,]$genes

gs<-kmes$genes[kmes$cols=="cyan"]
#cols<-cols[kmes$kme>=kME.threshold]
datExpr_cyan = datExpr[,gs]
head(datExpr_cyan[,1:8])
ncol(datExpr_cyan)
cyan_gene_list <- data.frame(colnames(datExpr_cyan))
#}#yellow module we have 3587 genes!!! 3951 genes for turquoise
###write into the file:
write.csv(cyan_gene_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_cyan_list_VST_distachyon.csv")

adj_mat<-adjacency(datExpr_cyan,power=9)
###
head(adj_mat[,1:8])
######### Get the top edges of the blue module:
top.n.edges = 3000
min.edge = 2
#head(adj_mat)
message("number of genes present = ", nrow(adj_mat))#number of genes present = 3587 Distachyon: 248
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)#adjacency.threshold = 0.690966673259787#distachyon:0.0952776666923627
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))#number of genes present = 158
#distachyon: 351
#adj_mat
network <- graph.adjacency(adj_mat,mode = "undirected",weighted = TRUE)
head(network)
E(network)
V(network)
deg <- degree(network, mode="all")
V(network)$vertex_degree <-  deg

network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
V(network)$color <- "cyan"
scale_factor <- 0.1
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#_cyan.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
#plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=2.5, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network,pch=19,arrow.mode=0,mark.border=NA, edge.color = "#d1d1e0", vertex.label="",vertex.size=V(network)$vertex_degree*scale_factor, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
dev.off()
