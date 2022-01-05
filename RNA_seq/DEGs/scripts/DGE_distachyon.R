#Written by Li Lei 12-19-2019 Berkeley
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("rlang")
biocLite("DESeq2")
biocLite("dplyr")
library(edgeR)
library(DESeq2)
library(dplyr)
require(dplyr)
#remove.packages("dplyr")
install.packages("dplyr")
#read the the count file for distachyon!!!
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Distachyon/New_counts_Feb2020/gene_counts/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#34310
#read the sorted meta data
timecourse <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/time_course_distachyon_sorted.txt",header = T)
head(timecourse)
nrow(timecourse)
#74
col.num <- which( colnames(counts) %in% timecourse$code)
#row.num <- which( metaData$code %in% colnames(counts) )
counts <- counts[, col.num]
ncol(counts)
#sort the data based on the rowname of the counts
counts <- counts %>% select(timecourse$code)
head(counts)
meta.data <- timecourse[,c(1,5,6,7)]
#mata <- metaData[row.num,]
colnames(counts) <- timecourse$sample_ID #change the headers
head(counts)

group <- paste0(meta.data$treatment, ".", meta.data$timepoint)
#meta.data
head(counts)
count.cutoff = 1 # 10
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]

nrow(counts)
#30717
ncol(counts)
#74
head(counts)
              
design <- model.matrix(~0 + group,meta.data)
#dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = ~ treatment+time+treatment:time)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = design)
nrow(meta.data)
ncol(counts)

ddsTC <- DESeq(dds)
result <- results(ddsTC)
head(result)44
#quality check
vsd <- vst(dds, blind=FALSE)
head(vsd)
colours <- c("#e6194b", "#3cb44b", "#ffe119", "#4363d8",
             "#f58231", "#911eb4", "#46f0f0", "#f032e6", 
             "#bcf60c", "#fabebe", "#008080", "#e6beff", 
             "#9a6324", "#fffac8", "#800000", "#aaffc3", 
             "#808000", "#ffd8b1", "#000075", "#808080")
#plotPCA(vsd,intgroup = c("treatment","timepoint"),col = colours)
pcaData <- plotPCA(vsd, intgroup=c("treatment","timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$timepoint <- as.factor(pcaData$timepoint)
ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=timepoint))+
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/quality_check_pca_dist.pdf", width = 15, height = 15)

#  ggplot(data = mpg) +
#    geom_point(mapping = aes(x = displ, y = hwy, color = cty)) 
head(vsd$timepoint)
head(vsd$treatment)
head(vsd$sizeFactor)
#heatmap:
library("pheatmap")
rld <- vst( ddsTC )
#head(heatmap.genes)
#expr.vst <- assay(DESeq2::vst(ddsTC))
head(rld)
library( "genefilter" )
#install.packages("gplots")
library(gplots)

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
head(topVarGenes)
#dev.off()
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heatmap_100_DGE_sylvaticum_vst.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='none', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", heat="#ff8000", drought="#ffff00",salt="#0066ff")[
             colData(rld)$treatment ] )
dev.off()

ir.pca <- prcomp(datExpr,
                 center = TRUE,
                 scale. = TRUE) 
print(ir.pca)
library(ggfortify)
#pdf(file = "/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_all_sample_pca_mad75.pdf",width = 10,height = 10)
autoplot(ir.pca,data = datTraits, colour = 'group',shape = FALSE, label.size = 3)
ggsave("/global/projectb/scratch/llei2019/CSP_Kranthi/VST_sus_all_sample_pca_mad75.pdf", width = 30, height = 25)
dev.off()


#Since the salt-24hs is werid, so I have to delete it and redo it!
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#36927

head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
counts <- subset(counts,select = c(1:60))#take the subset of the data
counts$Shoot_salt_24h.s3 <- NULL
#meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/coldata_abiotic_shoot_nosalt24hs3.txt",header = T)
head(meta.data)
group <- paste0(meta.data$treatment, ".", meta.data$time)

head(counts)
count.cutoff = 1 # 10
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]

nrow(counts)
ncol(counts)
head(counts)
#29910
#dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = as.data.frame(meta.data), design = ~group )
#dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = as.data.frame(meta.data), design = ~0 + group )
#ddsTC <- DESeq(dds)
#result <- results(ddsTC)
#resultsNames(ddsTC)
#res <- results(ddsTC, name="group_Shoot_drought_10h_vs_Shoot_CK_10h")               
#res <- results(ddsTC, name="group_Shoot_drought_24h_vs_Shoot_CK_24h")               

design <- model.matrix(~0 + group,meta.data)
#dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = ~ treatment+time+treatment:time)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = design)


#mm <- model.matrix(~time + treatment:time, meta.data)
#mm <- model.matrix(~treatment + time + time:treatment, meta.data)

#all.zero <- apply(mm, 2, function(x) all(x == 0))
#mm.full <- mm[, !all.zero]
#mm.reduced <- model.matrix(~time + treatment, meta.data)
#mm.reduced <- model.matrix(~treatment + time, meta.data)

#ddsTC <- DESeq(dds, full = mm.full, reduced = mm.reduced, test="LRT")
ddsTC <- DESeq(dds2)
result <- results(ddsTC)
head(result)
#heatmap:
library("pheatmap")
rld <- vst( ddsTC )
head(heatmap.genes)
#expr.vst <- assay(DESeq2::vst(ddsTC))
head(rld)
library( "genefilter" )
#install.packages("gplots")
library(gplots)

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
head(topVarGenes)
dev.off()
pdf(file = "/global/homes/l/llei2019/bscratch/RNAseq_Syl_sgordon/heatmap_100_DGE_nosalt24hs3.pdf",width = 15,height = 15)

par(mar=c(1,4,1,1))#trbl
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='none', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", heat="#ff8000", drought="#ffff00",salt="#0066ff")[
             colData(rld)$treatment ] )

#heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
#           trace="none", dendrogram='none', 
#           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
#           ColSideColors = c( CK="gray", heat="#ff8000", drought="#ffff00",salt="#0066ff")[
#             colData(rld)$treatment ] )
dev.off()

####check the outliers, Avinash suggested that we should use the lumi:

source("https://bioconductor.org/biocLite.R")
biocLite("lumi")
#detectOutlier(x, metric = "euclidean", standardize = TRUE, Th = 2, ifPlot = FALSE)

##control versus treatments
resultsNames(ddsTC)
#res_heat1h <- results(ddsTC, contrast=list("groupheat.1", "groupCK.1"))
res_heat1h <- na.omit(results(ddsTC, contrast=list("groupheat.1h", "groupCK.1h")))               
res_heat2h <- na.omit(results(ddsTC, contrast=list("groupheat.2h", "groupCK.2h")))               
res_heat5h <- na.omit(results(ddsTC, contrast=list("groupheat.5h", "groupCK.5h")))               
res_heat10h <- na.omit(results(ddsTC, contrast=list("groupheat.10h", "groupCK.10h")))               
res_heat24h <- na.omit(results(ddsTC, contrast=list("groupheat.24h", "groupCK.24h"))) 
#na.omit(res_heat1h)
sigres_heat1h_DGE <- res_heat1h[(res_heat1h$padj<0.05 & abs(res_heat1h$log2FoldChange) > 1),]
#heat1h
sigres_heat1h_up <- res_heat1h[(res_heat1h$padj<0.05 & res_heat1h$log2FoldChange>1),]
nb_heat1h_up <- nrow(sigres_heat1h_up)
#3325
sigres_heat1h_down <- res_heat1h[(res_heat1h$padj<0.05 & res_heat1h$log2FoldChange< -1),]
nb_heat1h_down <- nrow(sigres_heat1h_down)
#3354

#heat2h
sigres_heat2h_DGE <- res_heat2h[(res_heat2h$padj<0.05 & abs(res_heat2h$log2FoldChange)>1),]

sigres_heat2h_up <- res_heat2h[(res_heat2h$padj<0.05 & res_heat2h$log2FoldChange>1),]
nb_heat2h_up <- nrow(sigres_heat2h_up)
#2727
sigres_heat2h_down <- res_heat2h[(res_heat2h$padj<0.05 & res_heat2h$log2FoldChange< -1),]
nb_heat2h_down <- nrow(sigres_heat2h_down)
#3274
#heat5h
sigres_heat5h_DGE <- res_heat5h[(res_heat5h$padj<0.05 & abs(res_heat5h$log2FoldChange)>1),]
sigres_heat5h_up <- res_heat5h[(res_heat5h$padj<0.05 & res_heat5h$log2FoldChange>1),]
nb_heat5h_up <- nrow(sigres_heat5h_up)
#1912
sigres_heat5h_down <- res_heat5h[(res_heat5h$padj<0.05 & res_heat5h$log2FoldChange< -1),]
nb_heat5h_down <- nrow(sigres_heat5h_down)
#1866

#heat10h
sigres_heat10h_DGE <- res_heat10h[(res_heat10h$padj<0.05 & abs(res_heat10h$log2FoldChange)>1),]

sigres_heat10h_up <- res_heat10h[(res_heat10h$padj<0.05 & res_heat10h$log2FoldChange>1),]
nb_heat10h_up <- nrow(sigres_heat10h_up)
#2369
sigres_heat10h_down <- res_heat10h[(res_heat10h$padj<0.05 & res_heat10h$log2FoldChange< -1),]
nb_heat10h_down <- nrow(sigres_heat10h_down)
#2034

#heat24h
sigres_heat24h_up <- res_heat24h[(res_heat24h$padj<0.05 & abs(res_heat24h$log2FoldChange)>1),]
sigres_heat24h_up <- res_heat24h[(res_heat24h$padj<0.05 & res_heat24h$log2FoldChange>1),]
nb_heat24h_up <- nrow(sigres_heat24h_up)
#3612
sigres_heat24h_down <- res_heat24h[(res_heat24h$padj<0.05 & res_heat24h$log2FoldChange< -1),]
nb_heat24h_down <- nrow(sigres_heat24h_down)
#2586


####Drought
res_drought1h <- na.omit(results(ddsTC, contrast=list("groupdrought.1h", "groupCK.1h")))               
res_drought2h <- na.omit(results(ddsTC, contrast=list("groupdrought.2h", "groupCK.2h")))               
res_drought5h <- na.omit(results(ddsTC, contrast=list("groupdrought.5h", "groupCK.5h")))               
res_drought10h <- na.omit(results(ddsTC, contrast=list("groupdrought.10h", "groupCK.10h")))               
res_drought24h <- na.omit(results(ddsTC, contrast=list("groupdrought.24h", "groupCK.24h"))) 
#na.omit(res_heat1h)
#heat1h
sigres_drought1h_DGE <- res_drought1h[(res_drought1h$padj<0.05 & abs(res_drought1h$log2FoldChange)>1),]

sigres_drought1h_up <- res_drought1h[(res_drought1h$padj<0.05 & res_drought1h$log2FoldChange > 1),]
nb_drought1h_up <- nrow(sigres_drought1h_up)
#2479
sigres_drought1h_down <- res_drought1h[(res_drought1h$padj<0.05 & res_drought1h$log2FoldChange < -1),]
nb_drought1h_down <- nrow(sigres_drought1h_down)
#1397

#heat2h
sigres_drought2h_DGE <- res_drought2h[(res_drought2h$padj<0.05 & abs(res_drought2h$log2FoldChange)>1),]
sigres_drought2h_up <- res_drought2h[(res_drought2h$padj<0.05 & res_drought2h$log2FoldChange > 1),]
nb_drought2h_up <- nrow(sigres_drought2h_up)
#3652
sigres_drought2h_down <- res_drought2h[(res_drought2h$padj<0.05 & res_drought2h$log2FoldChange< -1),]
nb_drought2h_down <- nrow(sigres_drought2h_down)
#2675
#heat5h
sigres_drought5h_DGE <- res_drought5h[(res_drought5h$padj<0.05 & abs(res_drought5h$log2FoldChange)>1),]

sigres_drought5h_up <- res_drought5h[(res_drought5h$padj<0.05 & res_drought5h$log2FoldChange> 1),]
nb_drought5h_up <- nrow(sigres_drought5h_up)
#5434
sigres_drought5h_down <- res_drought5h[(res_drought5h$padj<0.05 & res_drought5h$log2FoldChange< -1),]
nb_drought5h_down <- nrow(sigres_drought5h_down)
#5355

#heat10h
sigres_drought10h_DGE <- res_drought10h[(res_drought10h$padj<0.05 & abs(res_drought10h$log2FoldChange)>1),]

sigres_drought10h_up <- res_drought10h[(res_drought10h$padj<0.05 & res_drought10h$log2FoldChange>1),]
nb_drought10h_up <- nrow(sigres_drought10h_up)
#5814
sigres_drought10h_down <- res_drought10h[(res_drought10h$padj<0.05 & res_drought10h$log2FoldChange< -1),]
nb_drought10h_down <- nrow(sigres_drought10h_down)
#5297

#heat24h
sigres_drought24h_DGE <- res_drought24h[(res_drought24h$padj<0.05 & abs(res_drought24h$log2FoldChange)>1),]

sigres_drought24h_up <- res_drought24h[(res_drought24h$padj<0.05 & res_drought24h$log2FoldChange>1),]
nb_drought24h_up <- nrow(sigres_drought24h_up)
#6341
sigres_drought24h_down <- res_drought24h[(res_drought24h$padj<0.05 & res_drought24h$log2FoldChange< -1),]
nb_drought24h_down <- nrow(sigres_drought24h_down)
#5921

####salt
res_salt1h <- na.omit(results(ddsTC, contrast=list("groupsalt.1h", "groupCK.1h")))               
res_salt2h <- na.omit(results(ddsTC, contrast=list("groupsalt.2h", "groupCK.2h")))               
res_salt5h <- na.omit(results(ddsTC, contrast=list("groupsalt.5h", "groupCK.5h")))               
res_salt10h <- na.omit(results(ddsTC, contrast=list("groupsalt.10h", "groupCK.10h")))               
res_salt24h <- na.omit(results(ddsTC, contrast=list("groupsalt.24h", "groupCK.24h"))) 
#na.omit(res_heat1h)
#heat1h
sigres_salt1h_DGE <- res_salt1h[(res_salt1h$padj<0.05 & abs(res_salt1h$log2FoldChange)>1),]

sigres_salt1h_up <- res_salt1h[(res_salt1h$padj<0.05 & res_salt1h$log2FoldChange > 1),]
nb_salt1h_up <- nrow(sigres_salt1h_up)
head(sigres_salt1h_up)
#1458
sigres_salt1h_down <- res_salt1h[(res_salt1h$padj<0.05 & res_salt1h$log2FoldChange< -1),]
nb_salt1h_down <- nrow(sigres_salt1h_down)
#186

#heat2h
sigres_salt2h_DGE <- res_salt2h[(res_salt2h$padj<0.05 & abs(res_salt2h$log2FoldChange)>1),]
sigres_salt2h_up <- res_salt2h[(res_salt2h$padj<0.05 & res_salt2h$log2FoldChange > 1),]
nb_salt2h_up <- nrow(sigres_salt2h_up)
#587
sigres_salt2h_down <- res_salt2h[(res_salt2h$padj<0.05 & res_salt2h$log2FoldChange < -1),]
nb_salt2h_down <- nrow(sigres_salt2h_down)
#134

#heat5h
sigres_salt5h_DGE <- res_salt5h[(res_salt5h$padj<0.05 & abs(res_salt5h$log2FoldChange)>1),]
sigres_salt5h_up <- res_salt5h[(res_salt5h$padj<0.05 & res_salt5h$log2FoldChange>1),]
nb_salt5h_up <- nrow(sigres_salt5h_up)
#296
sigres_salt5h_down <- res_salt5h[(res_salt5h$padj<0.05 & res_salt5h$log2FoldChange< -1),]
nb_salt5h_down <- nrow(sigres_salt5h_down)
#70

#heat10h
sigres_salt10h_DGE <- res_salt10h[(res_salt10h$padj<0.05 & abs(res_salt10h$log2FoldChange)>1),]

sigres_salt10h_up <- res_salt10h[(res_salt10h$padj<0.05 & res_salt10h$log2FoldChange>1),]
nb_salt10h_up <- nrow(sigres_salt10h_up)
#480
sigres_salt10h_down <- res_salt10h[(res_salt10h$padj<0.05 & res_salt10h$log2FoldChange< -1),]
nb_salt10h_down <- nrow(sigres_salt10h_down)
#288
#heat24h
sigres_salt24h_up <- res_salt24h[(res_salt24h$padj<0.05 & abs(res_salt24h$log2FoldChange)>1),]

sigres_salt24h_up <- res_salt24h[(res_salt24h$padj<0.05 & res_salt24h$log2FoldChange>1),]
nb_salt24h_up <- nrow(sigres_salt24h_up)
#789
sigres_salt24h_down <- res_salt24h[(res_salt24h$padj<0.05 & res_salt24h$log2FoldChange< -1),]
nb_salt24h_down <- nrow(sigres_salt24h_down)
#631



##GGPLOt
library(ggplot2)
heat <- data.frame(Categories=rep(c("Up", "Down"), each=5),
                   Time=rep(c("1h", "2h", "5h","10h","24h"),2),
                   Number=c(nb_heat1h_up, nb_heat2h_up, nb_heat5h_up, nb_heat10h_up, nb_heat24h_up, nb_heat1h_down, nb_heat2h_down, nb_heat5h_down, nb_heat10h_down, nb_heat24h_down))
head(heat)
str(heat)
heat$Time <- factor(heat$Time, levels = c("1h", "2h", "5h","10h","24h"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat_DEG_nb_nosalt24hs3_dist.pdf",width = 8,height = 8)
p <- ggplot(data=heat, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())

heatplot <- p + labs(title="Heat", 
                     x="Time (hours)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(heatplot)
dev.off()

#drought
drought <- data.frame(Categories=rep(c("Up", "Down"), each=5),
                      Time=rep(c("1h", "2h", "5h","10h","24h"),2),
                      Number=c(nb_drought1h_up, nb_drought2h_up, nb_drought5h_up, nb_drought10h_up, nb_drought24h_up, nb_drought1h_down, nb_drought2h_down, nb_drought5h_down, nb_drought10h_down, nb_drought24h_down))
head(drought)
str(drought)
drought$Time <- factor(drought$Time, levels = c("1h", "2h", "5h","10h","24h"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought_DEG_nb_nosalt_24hs3_dist.pdf",width = 8,height = 8)
p <- ggplot(data=drought, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())

droughtplot <- p + labs(title="Drought", 
                        x="Time (hours)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(droughtplot)
dev.off()

#salt
salt <- data.frame(Categories=rep(c("Up", "Down"), each=5),
                   Time=rep(c("1h", "2h", "5h","10h","24h"),2),
                   Number=c(nb_salt1h_up, nb_salt2h_up, nb_salt5h_up, nb_salt10h_up, nb_salt24h_up, nb_salt1h_down, nb_salt2h_down, nb_salt5h_down, nb_salt10h_down, nb_salt24h_down))
head(salt)
str(salt)
salt$Time <- factor(salt$Time, levels = c("1h", "2h", "5h","10h","24h"))
par(mfrow = c(1,1))
par(mar=c(5,5,5,5))

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt_DEG_nb_no_salt_24h_s3_dist.pdf",width = 8,height = 8)
p <- ggplot(data=salt, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())

saltplot <- p + labs(title="Salt", 
                     x="Time (hours)", y = "The number of genes")+
  scale_fill_manual(values=c('black','lightgray')) + theme_bw() +theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.text.x = element_text(color = "black", size = 20),
    axis.text.y = element_text(color = "black", size = 20),  
    axis.title.x = element_text(color = "black", size = 20),
    axis.title.y = element_text(color = "black", size = 20)
  )

print(saltplot)
dev.off()
##Take the union of the genes with up or downregulating 
#Take the union of the upregulating genes under heat tretment:
#upregulating:
heat1_2_up <- union(rownames(sigres_heat1h_up),rownames(sigres_heat2h_up))
heat1_5_up <- union(heat1_2_up,rownames(sigres_heat5h_up))
heat1_10_up <- union(heat1_5_up,rownames(sigres_heat10h_up))
heat1_24_up <- union(heat1_10_up,rownames(sigres_heat24h_up))
length(heat1_24_up)
#6540
#downregulating:
heat1_2_down <- union(rownames(sigres_heat1h_down),rownames(sigres_heat2h_down))
heat1_5_down <- union(heat1_2_down,rownames(sigres_heat5h_down))
heat1_10_down <- union(heat1_5_down,rownames(sigres_heat10h_down))
heat1_24_down <- union(heat1_10_down,rownames(sigres_heat24h_down))
length(heat1_24_down)
#6043
####Drought:
#Take the union of the upregulating genes under drought tretment:
#upregulating:
drought1_2_up <- union(rownames(sigres_drought1h_up),rownames(sigres_drought2h_up))
drought1_5_up <- union(drought1_2_up,rownames(sigres_drought5h_up))
drought1_10_up <- union(drought1_5_up,rownames(sigres_drought10h_up))
drought1_24_up <- union(drought1_10_up,rownames(sigres_drought24h_up))
length(drought1_24_up)
#8186
#downregulating:
drought1_2_down <- union(rownames(sigres_drought1h_down),rownames(sigres_drought2h_down))
drought1_5_down <- union(drought1_2_down,rownames(sigres_drought5h_down))
drought1_10_down <- union(drought1_5_down,rownames(sigres_drought10h_down))
drought1_24_down <- union(drought1_10_down,rownames(sigres_drought24h_down))
length(drought1_24_down)
#8353

#Salt:
#Take the union of the upregulating genes under salt tretment:
#upregulating:
salt1_2_up <- union(rownames(sigres_salt1h_up),rownames(sigres_salt2h_up))
salt1_5_up <- union(salt1_2_up,rownames(sigres_salt5h_up))
salt1_10_up <- union(salt1_5_up,rownames(sigres_salt10h_up))
salt1_24_up <- union(salt1_10_up,rownames(sigres_salt24h_up))
length(salt1_24_up)
#2503
#downregulating:
salt1_2_down <- union(rownames(sigres_salt1h_down),rownames(sigres_salt2h_down))
salt1_5_down <- union(salt1_2_down,rownames(sigres_salt5h_down))
salt1_10_down <- union(salt1_5_down,rownames(sigres_salt10h_down))
salt1_24_down <- union(salt1_10_down,rownames(sigres_salt24h_down))
length(salt1_24_down)
#1008

#Venn diagram:
#install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)
#myCol <- brewer.pal(3, "Purples")
myCol <- c("#ff8000","#ffff00","#0066ff")
#pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/Venn_up_no_salt_24h_s3.pdf",width = 8,height = 8)
#up
setwd("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/")
#dev.off()
venn.diagram(
  x = list(heat1_24_up, drought1_24_up, salt1_24_up),
  category.names = c("Heat" , "Drought " , "Salt"),
  filename = 'up-DGE-abiotic_dist.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 400,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
#down:
#up
venn.diagram(
  x = list(heat1_24_down, drought1_24_down, salt1_24_down),
  category.names = c("Heat" , "Drought " , "Salt"),
  filename = 'down-DGE-abiotic_dist.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 400,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

###Heatmap:


?brewer.pal
###Heatmap: 
heat1h <- results(ddsTC, contrast=list("groupheat.1h", "groupCK.1h"))              
heat2h <- results(ddsTC, contrast=list("groupheat.2h", "groupCK.2h"))               
heat5h <- results(ddsTC, contrast=list("groupheat.5h", "groupCK.5h"))               
heat10h <-results(ddsTC, contrast=list("groupheat.10h", "groupCK.10h"))               
heat24h <-results(ddsTC, contrast=list("groupheat.24h", "groupCK.24h"))

nrow(heat24h)
#29910
drought1h <- results(ddsTC, contrast=list("groupdrought.1h", "groupCK.1h"))              
drought2h <- results(ddsTC, contrast=list("groupdrought.2h", "groupCK.2h"))               
drought5h <- results(ddsTC, contrast=list("groupdrought.5h", "groupCK.5h"))               
drought10h <-results(ddsTC, contrast=list("groupdrought.10h", "groupCK.10h"))               
drought24h <-results(ddsTC, contrast=list("groupdrought.24h", "groupCK.24h"))

nrow(heat1h)#29910
salt1h <- results(ddsTC, contrast=list("groupsalt.1h", "groupCK.1h"))              
salt2h <- results(ddsTC, contrast=list("groupsalt.2h", "groupCK.2h"))               
salt5h <- results(ddsTC, contrast=list("groupsalt.5h", "groupCK.5h"))               
salt10h <-results(ddsTC, contrast=list("groupsalt.10h", "groupCK.10h"))               
salt24h <-results(ddsTC, contrast=list("groupsalt.24h", "groupCK.24h"))

combined <- data.frame(cbind (heat1h$log2FoldChange,heat1h$padj,
                              heat2h$log2FoldChange,heat2h$padj,
                              heat5h$log2FoldChange,heat5h$padj,
                              heat10h$log2FoldChange,heat10h$padj,
                              heat24h$log2FoldChange,heat24h$padj,
                              drought1h$log2FoldChange,drought1h$padj,
                              drought2h$log2FoldChange,drought2h$padj,
                              drought5h$log2FoldChange,drought5h$padj,
                              drought10h$log2FoldChange,drought10h$padj,
                              drought24h$log2FoldChange,drought24h$padj,
                              salt1h$log2FoldChange,salt1h$padj,
                              salt2h$log2FoldChange,salt2h$padj,
                              salt5h$log2FoldChange,salt5h$padj,
                              salt10h$log2FoldChange,salt10h$padj,
                              salt24h$log2FoldChange,salt24h$padj))
rownames(combined) <- rownames(salt24h)
colnames(combined) <- c("heat1h-log2FoldChange","heat1h-padj",
                        "heat2h-log2FoldChange","heat2h-padj",
                        "heat5h-log2FoldChange","heat5h-padj",
                        "heat10h-log2FoldChange","heat10h-padj",
                        "heat24h-log2FoldChange","heat24h-padj",
                        "drought1h-log2FoldChange","drought1h-padj",
                        "drought2h-log2FoldChange","drought2h-padj",
                        "drought5h-log2FoldChange","drought5h-padj",
                        "drought10h-log2FoldChange","drought10h-padj",
                        "drought24h-log2FoldChange","drought24h-padj",
                        "salt1h-log2FoldChange","salt1h-padj",
                        "salt2h-log2FoldChange","salt2h-padj",
                        "salt5h-log2FoldChange","salt5h-padj",
                        "salt10h-log2FoldChange","salt10h-padj",
                        "salt24h-log2FoldChange","salt24h-padj")
head(combined)

write.table(
  combined,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/DGE_shoot_abiotic_distachyon.txt",
  quote=F,
  row.names=T,
  sep="\t")



combined <- na.omit(combined)
nrow(combined)
DGE_all <- combined[(combined$`heat1h-padj` <0.05 &
                       combined$`heat2h-padj` <0.05 &
                       combined$`heat5h-padj` <0.05 &
                       combined$`heat10h-padj` <0.05 &
                       combined$`heat24h-padj` <0.05 &
                       combined$`drought1h-padj` <0.05 &
                       combined$`drought2h-padj` <0.05 &  
                       combined$`drought5h-padj` <0.05 &
                       combined$`drought10h-padj` <0.05 &
                       combined$`drought24h-padj` <0.05 &
                       combined$`salt1h-padj` <0.05 &
                       combined$`salt2h-padj` <0.05 &  
                       combined$`salt5h-padj` <0.05 &
                       combined$`salt10h-padj` <0.05 &
                       combined$`salt24h-padj` <0.05),]
nrow(DGE_all)
#7
head(DGE_all)
#?colData()
#res_1hdrought <- results(ddsTC, name="time1h.treatmentdrought")
#res <- results(ddsTC, contrast=list("treatmentdrought.time2h", "treatmentCK.time2h"))               
#res <- results(ddsTC, contrast=list("treatmentdrought", "treatmentCK"))   
#install.packages("openssl")
#install.packages("base64")
#install.packages("httr")
#install.packages("GEOquery")
#install.packages("biomaRt")
#source("https://bioconductor.org/biocLite.R")
#biocLite("lumi")
#biocLite("GenomicFeatures")
#library(lumi)
