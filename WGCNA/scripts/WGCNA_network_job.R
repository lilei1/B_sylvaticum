#!/usr/bin/env Rscript
#   Script to do the WGCNA analysis
# Written by Li Lei, Nov13, 2019 in Walnut Creek, CA.

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

install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
install.packages("igraph","ggfortify")
install.packages("ggfortify")
library(igraph)

#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute") 

library(textshape)
###Cleaning Data
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

tpm <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_normalized_count_#1.txt",sep="\t",header=T)
#use log(tpm+1) to do a little bit conversion
head(tpm)
str(tpm)
#tpm[1,]
#hist(tpm$Shoot_CK_1h.s1)
#femData <- sapply(tpm,function(x) log(x+1))
femData <- tpm
#column_to_rownames(tpm,"X")
#row.names(femData) <- tpm$GID
femData <-data.frame(femData)

### select genes with MAD >75%, at least MAD > 0.01 femData <- read.csv("/Users/LiLei/Downloads/tpm_counts.csv")
dim(femData)
head(femData)
m.mad <- apply(femData,1,mad)#median absolute deviation (MAD)
head(m.mad)
#femData <- femData[which(m.mad > 0),]
#This is to take the 75% genes
#femData <- femData[which(m.mad > 
 #                               max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
#above 85%
femData <- femData[which(m.mad > 
                               max(quantile(m.mad, probs=seq(0, 1, 0.05))[4],0.01)),]

###
datExpr0 <- as.data.frame(t(femData))
head(datExpr0[,1:8])
dim(datExpr0)
#[1]    60 21663
#[1]    60 24742
#[1]    60 24552 85%
rownames(datExpr0)
colnames(datExpr0)
#names(datExpr0) = femData$G_id
#rownames(datExpr0) = names(femData)[-c(1)]
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
nGenes = ncol(datExpr0)#24552
nSamples = nrow(datExpr0)
head(datExpr0[,1:8])
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small. sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_abiotic4#.pdf",width = 25,height = 20)
plot(sampleTree,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#do a PCA analysis for the samples to see what will happen:
ir.pca <- prcomp(datExpr0,
                 center = TRUE,
                 scale. = TRUE) 
print(ir.pca)
meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
head(meta.data)
head(ir.pca)
head(datExpr0)

library(ggfortify)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_abiotic_forWGCNA_pca.pdf",width = 10,height = 10)
autoplot(ir.pca,data = meta.data, colour = 'group')
dev.off()
#####
###Automatic block-wise network
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
head(sft)
power = sft$powerEstimate
power
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/power_soft_thresh_VST_abiotic#3.pdf",width = 10,height = 10)
# Scale-free topology fit index as a function of the soft-thresholding power plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = "Scale independence")

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")
abline(h=0.8,col="red")
# this line corresponds to using an R^2 cut-off of h abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#abline(h=0.85,col="red")
dev.off()

#power=8 seems good!
bwnet = blockwiseModules(datExpr0,maxBlockSize = nGenes,
                         power = 8, TOMType = type, minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = "allAin_TOM-blockwise", verbose = 3)
head(bwnet)
bwLabels = bwnet$colors
table(bwnet$colors)
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
#558 8923 3587 2721 2579 2342  857  834  359  350  291  250  201  174  138  136  118   88   46 
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)


#?matchLabels

# open a graphics window
sizeGrWindow(6,6)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/Dendro_abiotic_VST.pdf",width = 25,height = 25)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

head(bwnet)
MEs = bwnet$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/eigennetwork_abiotic_VST.pdf",width = 15,height = 13)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
#MM = as.data.frame(cor(datExpr0, MEs, use = "p"))
#kIM = intramodularConnectivity(TOM, bwModuleColors, scaleByMax = TRUE) 
#head(kIM)
#tophub <- chooseTopHubInEachModule(datExpr0,bwModuleColors)
#head(tophub)
#?chooseTopHubInEachModule
#?signedKME
#extract top 50 hub genes for each module!!!
MEs <- moduleEigengenes(datExpr0, bwModuleColors)$eigengenes
head(MEs)
kMEs <- signedKME(datExpr0, MEs)

# rank the genes for each module on kMEs
rankGenes <- function(x){
  kMErank <- rank(-kMEs[ ,x])
  genes <- rownames(kMEs)
  genes <- genes[order(kMErank)]
  genes[1:50]#top 50 hub genes!!!
}

topGenes <- lapply(1:ncol(kMEs), rankGenes)

# Get the top results in a data.frame
topGenes <- do.call(cbind, topGenes)
colnames(topGenes) <- substr(colnames(kMEs), start=4, stop=30)
hubgenes <- data.frame(topGenes)
write.csv(hubgenes, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_50hubgenes_list.csv")

######### Get the top edges:
top.n.edges = 3000
min.edge = 2
adj_mat<-adjacency(datExpr0,power=8)
head(adj_mat)
message("number of genes present = ", nrow(adj_mat))
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))
head(adj_mat)
network <- graph.adjacency(adj_mat)
#head(network)
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
plot(network, vertex.color="yellow",vertex.label="",vertex.size=2, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

dev.off()
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
##It works! But how can I figure out the modules?
############Subset the network:
top.n.edges = 3000
min.edge = 2
adj_mat<-adjacency(datExpr0,power=8)

probes = colnames(datExpr0)
modules = c("yellow")
inModule = is.finite(match(bwModuleColors, modules))
modProbes = probes[inModule]
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];

modTOM = adj_mat[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
#plot the yellow module!!!
message("number of genes present = ", nrow(modTOM))
if(!is.na(top.n.edges) & nrow(adj_mat) >= 100){
  adjacency.threshold = sort(as.numeric(adj_mat), decreasing=T)[(top.n.edges*2)+nrow(adj_mat)]
}
if(nrow(adj_mat) < 100)
  adjacency.threshold = 0.01

message("adjacency.threshold = ", adjacency.threshold)
adj_mat[adj_mat > adjacency.threshold] <- 1
adj_mat[adj_mat < adjacency.threshold] <- 0
diag(adj_mat) <- 0
rs<-rowSums(adj_mat)
if(verbose) cat("removing unconnected nodes\n")
adj_mat<-adj_mat[rs>min.edge,rs>min.edge]
message("number of genes present = ", nrow(adj_mat))
head(adj_mat)
network <- graph.adjacency(adj_mat)
#head(network)
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_all_abiotic_VST#.pdf",width = 25,height = 25)
#plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)
plot(network, vertex.color="yellow",vertex.label="",vertex.size=2, layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)

dev.off()


######It works!!!
#####Build the network!!!
adj <- modTOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network) 
#######!!!!!!!!!!!!!!!!!!



###Tomplot
#load(bwnet$TOMFiles[1], verbose=T)
#TOM <- as.matrix(TOM)
#dissTOM = 1-TOM
#plotTOM = dissTOM^7
#diag(plotTOM) = NA
#pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/Tomplot_abiotic_voom.pdf",width = 25,height = 25)
#TOMplot(plotTOM, bwnet$dendrograms, bwModuleColors, 
#        main = "Network heatmap plot, all genes")
#dev.off()

###subset the gens and test partial of the data:
load(bwnet$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
nSelect = 400
nGenes = ncol(datExpr0)
set.seed(10);
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = bwModuleColors[select]
#sizeGrWindow(9,9)
plotDiss = selectTOM^7
diag(plotDiss) = NA
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/Tomplot_abiotic_400gene_VST.pdf",width = 25,height = 25)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()
#######It works!!!!

###network
#probes = colnames(datExpr0)
#dimnames(TOM) <- list(probes, probes)

#cyt = exportNetworkToCytoscape(TOM,
#                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
#                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
#                               weighted = TRUE, threshold = 0,
#                               nodeNames = probes, nodeAttr = bwModuleColors)

###Subset the network:
probes = colnames(datExpr0)
modules = c("yellow")
inModule = is.finite(match(bwModuleColors, modules))
modProbes = probes[inModule]
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
#head(TOM)
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("abiotic-yellow-edges-", paste(modules, collapse="-"), ".txt", sep=""),
nodeFile = paste("abiotic-yellow-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.02,
nodeNames = modProbes,
#altNodeNames = modGenes,
nodeAttr = bwModuleColors[inModule])
getwd()
######It works!!!
#####Build the network!!!
adj <- modTOM
adj[adj > 0.1] = 1
adj[adj != 1] = 0
network <- graph.adjacency(adj)
network <- simplify(network)  # removes self-loops
#results <- blockwiseModules(data, power=6, TOMType="unsigned", networkType="unsigned")
#V(network)$color <- results$colors
par(mar=c(0,0,0,0))
# remove unconnected nodes
network <- delete.vertices(network, degree(network)==0)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/network_yellow_abiotic_voom.pdf",width = 25,height = 25)
plot(network, vertex.color="yellow", layout=layout.fruchterman.reingold(network), edge.arrow.size = 0.2)
plot(network, vertex.color="yellow", vertex.label="",vertex.size=2,layout=layout_in_circle(network), edge.arrow.size = 0.2)


dev.off()

#######It works!!!



#associate trait
trait <- "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_stress.txt"
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(datExpr0)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
robustY = ifelse(corType=="pearson",T,F)
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

sizeGrWindow(9, 5)
par(mar = c(5, 20, 5, 5),mfrow = c(1,1))
cex1 = 1
pdf("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pheno-module_abiotic_VST.pdf",width = 15,height = 15)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(datExpr0, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExpr0, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(datExpr0, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExpr0, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


module = "yellow"
pheno = "heat"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = bwModuleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
#??geneTraitSignificance
pdf("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/MEyellow_heat_VST.pdf",width = 15,height = 15)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
# Plot the dendrogram and the module colors underneath for block 2
#plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
#                    "Module colors", main = "Gene dendrogram and module colors in block 2",
#                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Plot the dendrogram and the module colors underneath for block 3
#plotDendroAndColors(bwnet$dendrograms[[3]], bwModuleColors[bwnet$blockGenes[[3]]],
#                    "Module colors", main = "Gene dendrogram and module colors in block 3",
#                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

#sizeGrWindow(12,9)
#plotDendroAndColors(geneTree,
#                    cbind(moduleColors, bwModuleColors),
#                    c("Single block", "2 blocks"),
#                    main = "Single block gene dendrogram and module colors", dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)


#cyan_list <- names(datExpr0)[bwModuleColors=="cyan"]
#write.csv(cyan_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_cyan_list_VST.csv")

#lightcyan_list <- names(datExpr0)[bwModuleColors=="lightcyan"]
#write.csv(lightcyan_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_lightcyan_list_VST.csv")

#red_list <- names(datExpr0)[bwModuleColors=="red"]
#write.csv(red_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_red_list_VST#.csv")

#black_list <- names(datExpr0)[bwModuleColors=="black"]
#write.csv(black_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_black_list_VST#.csv")

#magenta_list <- names(datExpr0)[bwModuleColors=="magenta"]
#write.csv(magenta_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_magenta_list_VST.csv")

#blue_list <- names(datExpr0)[bwModuleColors=="blue"]
#write.csv(blue_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_blue_list_VST.csv")

#midnightblue_list <- names(datExpr0)[bwModuleColors=="midnightblue"]
#write.csv(midnightblue_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_midnightblue_list_VST.csv")

#pink_list <- names(datExpr0)[bwModuleColors=="pink"]
#write.csv(pink_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_pink_list_VST#.csv")

#brown_list <- names(datExpr0)[bwModuleColors=="brown"]
#write.csv(brown_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_brown_list_VST.csv")

#green_list <- names(datExpr0)[bwModuleColors=="green"]
#write.csv(brown_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_green_list_VST.csv")

#greenyellow_list <- names(datExpr0)[bwModuleColors=="greenyellow"]
#write.csv(greenyellow_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_greenyellow_list_VST.csv")

#lightgreen_list <- names(datExpr0)[bwModuleColors=="lightgreen"]
#write.csv(lightgreen_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_lightgreen_list_VST.csv")

#yellow_list <- names(datExpr0)[bwModuleColors=="yellow"]
#write.csv(yellow_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_yellow_list_VST#.csv")

#grey_list <- names(datExpr0)[bwModuleColors=="grey"]
#write.csv(grey_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_grey_list_VST.csv")

#grey60_list <- names(datExpr0)[bwModuleColors=="grey60"]
#write.csv(grey60_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_grey60_list_VST.csv")


#turquoise_list <- names(datExpr0)[bwModuleColors=="turquoise"]
#write.csv(turquoise_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_turquoise_list_VST#.csv")

#salmon_list <- names(datExpr0)[bwModuleColors=="salmon"]
#write.csv(black_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_salmon_list_VST#.csv")

#purple_list <- names(datExpr0)[bwModuleColors=="purple"]
#write.csv(black_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_purple_list_VST#.csv")

#tan_list <- names(datExpr0)[bwModuleColors=="tan"]
#write.csv(black_list, file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_tan_list_VST#.csv")
#About the visualization of the network with igraph
