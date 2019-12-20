#Written by Li Lei 12-19-2019 Berkeley
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("rlang")
biocLite("DESeq2")
library(edgeR)
library(DESeq2)
#packageVersion("rlang")
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#36927

head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
counts <- subset(counts,select = c(1:60))#take the subset of the data

#meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/coldata_abiotic_shoot.txt",header = T)
head(meta.data)
group <- paste0(meta.data$treatment, ".", meta.data$time)

head(counts)
count.cutoff = 1 # 10
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]

nrow(counts)
#29952
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

head(dds)
#mm <- model.matrix(~time + treatment:time, meta.data)
#mm <- model.matrix(~treatment + time + time:treatment, meta.data)

#all.zero <- apply(mm, 2, function(x) all(x == 0))
#mm.full <- mm[, !all.zero]
#mm.reduced <- model.matrix(~time + treatment, meta.data)
#mm.reduced <- model.matrix(~treatment + time, meta.data)

#ddsTC <- DESeq(dds, full = mm.full, reduced = mm.reduced, test="LRT")
ddsTC <- DESeq(dds)
result <- results(ddsTC)
head(result)
#heatmap:
library("pheatmap")
rld <- rlog( ddsTC )
head(heatmap.genes)
#expr.vst <- assay(DESeq2::vst(ddsTC))
head(rld)
library( "genefilter" )
install.packages("gplots")
library(gplots)

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
head(topVarGenes)
dev.off()
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heatmap_100_DGE#.pdf",width = 15,height = 15)

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

#heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
#           trace="none", dendrogram='none', 
#           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
#           ColSideColors = c( CK="gray", heat="#ff8000", drought="#ffff00",salt="#0066ff")[
#             colData(rld)$treatment ] )
dev.off()
##control versus treatments
resultsNames(ddsTC)
res_heat1h <- na.omit(results(ddsTC, contrast=list("groupheat.1h", "groupCK.1h")))               
res_heat2h <- na.omit(results(ddsTC, contrast=list("groupheat.2h", "groupCK.2h")))               
res_heat5h <- na.omit(results(ddsTC, contrast=list("groupheat.5h", "groupCK.5h")))               
res_heat10h <- na.omit(results(ddsTC, contrast=list("groupheat.10h", "groupCK.10h")))               
res_heat24h <- na.omit(results(ddsTC, contrast=list("groupheat.24h", "groupCK.24h"))) 
#na.omit(res_heat1h)
#heat1h
sigres_heat1h_up <- res_heat1h[(res_heat1h$padj<0.05 & res_heat1h$log2FoldChange>0),]
nb_heat1h_up <- nrow(sigres_heat1h_up)
#5147
sigres_heat1h_down <- res_heat1h[(res_heat1h$padj<0.05 & res_heat1h$log2FoldChange<0),]
nb_heat1h_down <- nrow(sigres_heat1h_down)
#5047

#heat2h
sigres_heat2h_up <- res_heat2h[(res_heat2h$padj<0.05 & res_heat2h$log2FoldChange>0),]
nb_heat2h_up <- nrow(sigres_heat2h_up)
#5028
sigres_heat2h_down <- res_heat2h[(res_heat1h$padj<0.05 & res_heat2h$log2FoldChange<0),]
nb_heat2h_down <- nrow(sigres_heat2h_down)
#4878
#heat5h
sigres_heat5h_up <- res_heat5h[(res_heat5h$padj<0.05 & res_heat5h$log2FoldChange>0),]
nb_heat5h_up <- nrow(sigres_heat5h_up)
#3866
sigres_heat5h_down <- res_heat5h[(res_heat5h$padj<0.05 & res_heat5h$log2FoldChange<0),]
nb_heat5h_down <- nrow(sigres_heat5h_down)
#3362
#heat10h
sigres_heat10h_up <- res_heat10h[(res_heat10h$padj<0.05 & res_heat10h$log2FoldChange>0),]
nb_heat10h_up <- nrow(sigres_heat10h_up)
#4285
sigres_heat10h_down <- res_heat10h[(res_heat10h$padj<0.05 & res_heat10h$log2FoldChange<0),]
nb_heat10h_down <- nrow(sigres_heat10h_down)
#3662
#heat24h
sigres_heat24h_up <- res_heat24h[(res_heat24h$padj<0.05 & res_heat24h$log2FoldChange>0),]
nb_heat24h_up <- nrow(sigres_heat24h_up)
#5653
sigres_heat24h_down <- res_heat24h[(res_heat24h$padj<0.05 & res_heat24h$log2FoldChange<0),]
nb_heat24h_down <- nrow(sigres_heat24h_down)
#4653


####Drought
res_drought1h <- na.omit(results(ddsTC, contrast=list("groupdrought.1h", "groupCK.1h")))               
res_drought2h <- na.omit(results(ddsTC, contrast=list("groupdrought.2h", "groupCK.2h")))               
res_drought5h <- na.omit(results(ddsTC, contrast=list("groupdrought.5h", "groupCK.5h")))               
res_drought10h <- na.omit(results(ddsTC, contrast=list("groupdrought.10h", "groupCK.10h")))               
res_drought24h <- na.omit(results(ddsTC, contrast=list("groupdrought.24h", "groupCK.24h"))) 
#na.omit(res_heat1h)
#heat1h
sigres_drought1h_up <- res_drought1h[(res_drought1h$padj<0.05 & res_drought1h$log2FoldChange>0),]
nb_drought1h_up <- nrow(sigres_drought1h_up)
#3140
sigres_drought1h_down <- res_drought1h[(res_drought1h$padj<0.05 & res_drought1h$log2FoldChange<0),]
nb_drought1h_down <- nrow(sigres_drought1h_down)
#2068

#heat2h
sigres_drought2h_up <- res_drought2h[(res_drought2h$padj<0.05 & res_drought2h$log2FoldChange>0),]
nb_drought2h_up <- nrow(sigres_drought2h_up)
#5041
sigres_drought2h_down <- res_drought2h[(res_drought2h$padj<0.05 & res_drought2h$log2FoldChange<0),]
nb_drought2h_down <- nrow(sigres_drought2h_down)
#4367
#heat5h
sigres_drought5h_up <- res_drought5h[(res_drought5h$padj<0.05 & res_drought5h$log2FoldChange>0),]
nb_drought5h_up <- nrow(sigres_drought5h_up)
#7578
sigres_drought5h_down <- res_drought5h[(res_drought5h$padj<0.05 & res_drought5h$log2FoldChange<0),]
nb_drought5h_down <- nrow(sigres_drought5h_down)
#7462
#heat10h
sigres_drought10h_up <- res_drought10h[(res_drought10h$padj<0.05 & res_drought10h$log2FoldChange>0),]
nb_drought10h_up <- nrow(sigres_drought10h_up)
#7810
sigres_drought10h_down <- res_drought10h[(res_drought10h$padj<0.05 & res_drought10h$log2FoldChange<0),]
nb_drought10h_down <- nrow(sigres_drought10h_down)
#7421
#heat24h
sigres_drought24h_up <- res_drought24h[(res_drought24h$padj<0.05 & res_drought24h$log2FoldChange>0),]
nb_drought24h_up <- nrow(sigres_drought24h_up)
#8171
sigres_drought24h_down <- res_drought24h[(res_drought24h$padj<0.05 & res_drought24h$log2FoldChange<0),]
nb_drought24h_down <- nrow(sigres_drought24h_down)
#7770

####salt
res_salt1h <- na.omit(results(ddsTC, contrast=list("groupsalt.1h", "groupCK.1h")))               
res_salt2h <- na.omit(results(ddsTC, contrast=list("groupsalt.2h", "groupCK.2h")))               
res_salt5h <- na.omit(results(ddsTC, contrast=list("groupsalt.5h", "groupCK.5h")))               
res_salt10h <- na.omit(results(ddsTC, contrast=list("groupsalt.10h", "groupCK.10h")))               
res_salt24h <- na.omit(results(ddsTC, contrast=list("groupsalt.24h", "groupCK.24h"))) 
#na.omit(res_heat1h)
#heat1h
sigres_salt1h_up <- res_salt1h[(res_salt1h$padj<0.05 & res_salt1h$log2FoldChange>0),]
nb_salt1h_up <- nrow(sigres_salt1h_up)
head(sigres_salt1h_up)
#1868
sigres_salt1h_down <- res_salt1h[(res_salt1h$padj<0.05 & res_salt1h$log2FoldChange<0),]
nb_salt1h_down <- nrow(sigres_salt1h_down)
#424

#heat2h
sigres_salt2h_up <- res_salt2h[(res_salt2h$padj<0.05 & res_salt2h$log2FoldChange>0),]
nb_salt2h_up <- nrow(sigres_salt2h_up)
#679
sigres_salt2h_down <- res_salt2h[(res_salt2h$padj<0.05 & res_salt2h$log2FoldChange<0),]
nb_salt2h_down <- nrow(sigres_salt2h_down)
#297
#heat5h
sigres_salt5h_up <- res_salt5h[(res_salt5h$padj<0.05 & res_salt5h$log2FoldChange>0),]
nb_salt5h_up <- nrow(sigres_salt5h_up)
#294
sigres_salt5h_down <- res_salt5h[(res_salt5h$padj<0.05 & res_salt5h$log2FoldChange<0),]
nb_salt5h_down <- nrow(sigres_salt5h_down)
#211
#heat10h
sigres_salt10h_up <- res_salt10h[(res_salt10h$padj<0.05 & res_salt10h$log2FoldChange>0),]
nb_salt10h_up <- nrow(sigres_salt10h_up)
#668
sigres_salt10h_down <- res_salt10h[(res_salt10h$padj<0.05 & res_salt10h$log2FoldChange<0),]
nb_salt10h_down <- nrow(sigres_salt10h_down)
#435
#heat24h
sigres_salt24h_up <- res_salt24h[(res_salt24h$padj<0.05 & res_salt24h$log2FoldChange>0),]
nb_salt24h_up <- nrow(sigres_salt24h_up)
#4967
sigres_salt24h_down <- res_salt24h[(res_salt24h$padj<0.05 & res_salt24h$log2FoldChange<0),]
nb_salt24h_down <- nrow(sigres_salt24h_down)
#2925



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

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat_DEG_nb.pdf",width = 8,height = 8)
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

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought_DEG_nb.pdf",width = 8,height = 8)
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

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt_DEG_nb.pdf",width = 8,height = 8)
p <- ggplot(data=salt, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())

saltplot <- p + labs(title="Drought", 
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
#10002
#downregulating:
heat1_2_down <- union(rownames(sigres_heat1h_down),rownames(sigres_heat2h_down))
heat1_5_down <- union(heat1_2_down,rownames(sigres_heat5h_down))
heat1_10_down <- union(heat1_5_down,rownames(sigres_heat10h_down))
heat1_24_down <- union(heat1_10_down,rownames(sigres_heat24h_down))
length(heat1_24_down)
#8659
####Drought:
#Take the union of the upregulating genes under drought tretment:
#upregulating:
drought1_2_up <- union(rownames(sigres_drought1h_up),rownames(sigres_drought2h_up))
drought1_5_up <- union(drought1_2_up,rownames(sigres_drought5h_up))
drought1_10_up <- union(drought1_5_up,rownames(sigres_drought10h_up))
drought1_24_up <- union(drought1_10_up,rownames(sigres_drought24h_up))
length(drought1_24_up)
#10689
#downregulating:
drought1_2_down <- union(rownames(sigres_drought1h_down),rownames(sigres_drought2h_down))
drought1_5_down <- union(drought1_2_down,rownames(sigres_drought5h_down))
drought1_10_down <- union(drought1_5_down,rownames(sigres_drought10h_down))
drought1_24_down <- union(drought1_10_down,rownames(sigres_drought24h_down))
length(drought1_24_down)
#11041

#Salt:
#Take the union of the upregulating genes under salt tretment:
#upregulating:
salt1_2_up <- union(rownames(sigres_salt1h_up),rownames(sigres_salt2h_up))
salt1_5_up <- union(salt1_2_up,rownames(sigres_salt5h_up))
salt1_10_up <- union(salt1_5_up,rownames(sigres_salt10h_up))
salt1_24_up <- union(salt1_10_up,rownames(sigres_salt24h_up))
length(salt1_24_up)
#6475
#downregulating:
salt1_2_down <- union(rownames(sigres_salt1h_down),rownames(sigres_salt2h_down))
salt1_5_down <- union(salt1_2_down,rownames(sigres_salt5h_down))
salt1_10_down <- union(salt1_5_down,rownames(sigres_salt10h_down))
salt1_24_down <- union(salt1_10_down,rownames(sigres_salt24h_down))
length(salt1_24_down)
#3527

#Venn diagram:
install.packages("VennDiagram")
library(VennDiagram)
library(RColorBrewer)
#myCol <- brewer.pal(3, "Purples")
myCol <- c("#ff8000","#ffff00","#0066ff")
#up
venn.diagram(
  x = list(heat1_24_up, drought1_24_up, salt1_24_up),
  category.names = c("Heat" , "Drought " , "Salt"),
  filename = 'up-DGE-abiotic.png',
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
  filename = 'down-DGE-abiotic.png',
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
drought1h <- results(ddsTC, contrast=list("groupdrought.1h", "groupCK.1h"))              
drought2h <- results(ddsTC, contrast=list("groupdrought.2h", "groupCK.2h"))               
drought5h <- results(ddsTC, contrast=list("groupdrought.5h", "groupCK.5h"))               
drought10h <-results(ddsTC, contrast=list("groupdrought.10h", "groupCK.10h"))               
drought24h <-results(ddsTC, contrast=list("groupdrought.24h", "groupCK.24h"))

nrow(heat1h)
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
?colData()
#res_1hdrought <- results(ddsTC, name="time1h.treatmentdrought")
#res <- results(ddsTC, contrast=list("treatmentdrought.time2h", "treatmentCK.time2h"))               
#res <- results(ddsTC, contrast=list("treatmentdrought", "treatmentCK"))               
