#Written by Li Lei 12-19-2019 Berkeley
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#biocLite("httr", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
biocLite("Hmisc")
biocLite("DESeq2")
#biocLite("devtools")
#devtools::install_github('kevinblighe/EnhancedVolcano')

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager', dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#BiocManager::install('DESeq2', dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#BiocManager::install('latticeExtra', dependencies=TRUE, INSTALL_opts = c('--no-lock'))


#BiocManager::install('EnhancedVolcano', dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#BiocManager::install('ggplot2', dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#library(EnhancedVolcano)
#library(airway)
#library(magrittr)
library(ggplot2)
library(edgeR)
library(DESeq2)
library(tidyr)
library(ggrepel)
library(EnhancedVolcano)

#packageVersion("rlang")
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#36927

head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
#Just focus on the shoot along the time point
counts <- subset(counts,select = c(1:60))#take the subset of the data
#counts$Shoot_salt_24h.s3 <- NULL
#meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/coldata_abiotic_shoot.txt",header = T)
head(meta.data)
group <- paste0(meta.data$treatment, ".", meta.data$time)
#Find the best cutoff for the downstream analysis
#This is to convert a dataframe into single column.
counts_long <- gather(counts,code,counts)
head(counts_long)
#nrow(counts_long)
#this is to find the reasonable cutoff for the threshold setting
quantile(counts_long$counts,probs=seq(0,1,0.01))
#  0%     1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11%    12% 
#0      0      0      0      0      0      0      0      0      0      0      0      0 
#13%    14%    15%    16%    17%    18%    19%    20%    21%    22%    23%    24%    25% 
#0      0      0      0      0      0      0      0      0      0      0      0      0 
#26%    27%    28%    29%    30%    31%    32%    33%    34%    35%    36%    37%    38% 
#0      0      0      0      0      0      0      1      1      1      1      1      2 
#39%    40%    41%    42%    43%    44%    45%    46%    47%    48%    49%    50%    51% 
#2      3      3      4      4      5      6      7      8     10     11     13     15 
#52%    53%    54%    55%    56%    57%    58%    59%    60%    61%    62%    63%    64% 
#17     20     23     26     29     33     38     42     48     53     60     67     74 
#65%    66%    67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
#82     91    100    111    122    134    147    160    175    191    209    228    248 
#78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88%    89%    90% 
#270    295    320    349    380    414    452    494    540    592    651    717    794 
#91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
#  885    994   1129   1301   1534   1862   2375   3322   5744 569218 
head(counts)
count.cutoff = 5 # 56% of the data were above 5
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#[1] 25072
ncol(counts)
head(counts)

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
rld <- vst( ddsTC )
#head(heatmap.genes)
#expr.vst <- assay(DESeq2::vst(ddsTC))
head(rld)
library( "genefilter" )
#install.packages("gtools", "gdata", "caTools")
library(gplots)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
head(topVarGenes)
dev.off()
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heatmap_100_DGE_vst_syl.pdf",width = 15,height = 15)

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
####do filtering with raw count
count.cutoff = 5 # 10 #above56% of the datapoint is greater than 5;
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#25039
#CPM filtering
normalized.counts <- cpm(counts)
head(normalized.counts)
###below is for find the propriate threshold to do filtering based on the cpm normalization
normalized.counts.re <- data.frame(normalized.counts)
head(normalized.counts.re)
normalized.counts_long <- gather(normalized.counts.re,code,counts)
head(normalized.counts_long)
#nrow(counts_long)
#this is to find the reasonable cutoff for the threshold setting
quantile(normalized.counts_long$counts,probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 4.995701e-02 
#7%           8%           9%          10%          11%          12%          13% 
#5.911268e-02 6.783174e-02 8.479778e-02 1.115156e-01 1.343963e-01 1.573204e-01 1.814071e-01 
#14%          15%          16%          17%          18%          19%          20% 
#2.088814e-01 2.414296e-01 2.727203e-01 3.078148e-01 3.417429e-01 3.835590e-01 4.296072e-01 
#21%          22%          23%          24%          25%          26%          27% 
#4.823804e-01 5.322693e-01 5.927993e-01 6.538771e-01 7.235706e-01 7.941087e-01 8.812741e-01 
#28%          29%          30%          31%          32%          33%          34% 
#9.647608e-01 1.059398e+00 1.161570e+00 1.273716e+00 1.394981e+00 1.529318e+00 1.668347e+00 
#35%          36%          37%          38%          39%          40%          41% 
#1.814958e+00 1.979327e+00 2.150342e+00 2.331505e+00 2.535011e+00 2.739851e+00 2.968990e+00 
#42%          43%          44%          45%          46%          47%          48% 
#3.210356e+00 3.462148e+00 3.734636e+00 4.024631e+00 4.321101e+00 4.647758e+00 4.984598e+00 
#49%          50%          51%          52%          53%          54%          55% 
#5.342260e+00 5.711845e+00 6.117675e+00 6.534358e+00 6.974907e+00 7.445696e+00 7.931625e+00 
#56%          57%          58%          59%          60%          61%          62% 
#8.439328e+00 8.972948e+00 9.550031e+00 1.014911e+01 1.077719e+01 1.144304e+01 1.214145e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.288603e+01 1.365663e+01 1.447113e+01 1.532997e+01 1.623603e+01 1.718678e+01 1.819329e+01 
#70%          71%          72%          73%          74%          75%          76% 
#1.926496e+01 2.039283e+01 2.158683e+01 2.287661e+01 2.424809e+01 2.570640e+01 2.723837e+01 
#77%          78%          79%          80%          81%          82%          83% 
#2.889371e+01 3.067807e+01 3.256865e+01 3.462585e+01 3.680749e+01 3.917049e+01 4.175782e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.464346e+01 4.780005e+01 5.129312e+01 5.523629e+01 5.972040e+01 6.494892e+01 7.097583e+01 
#91%          92%          93%          94%          95%          96%          97% 
#7.815822e+01 8.696257e+01 9.770888e+01 1.115630e+02 1.298731e+02 1.567913e+02 1.991155e+02 
#98%          99%         100% 
#2.738284e+02 4.741552e+02 2.980236e+04  
count.cutoff = 1 #90% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#21794
nrow(counts)
#25039
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#21794
#counts[21822,]
grep("Brasy3G282900", rownames(counts))#check if this gene got filtered
#8339
#it suggested that the filtering seems good!
#21794/36927=0.5901915 only keep around 60% genes

## VST instead of voom
#counts[is.na(counts)] <- 0
#creat a group 
design <- model.matrix(~0 + group,meta.data)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = design)
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
#str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#This file is for WGCNA
write.table(x = expr.vst,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/normalized_VST_final_shoot_abiotic_formal.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
            )

###Do DEG
#design <- model.matrix(~0 + group,meta.data)
#dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = ~ treatment+time+treatment:time)
#dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta.data, design = design)

#head(dds)
#mm <- model.matrix(~time + treatment:time, meta.data)
#mm <- model.matrix(~treatment + time + time:treatment, meta.data)

#ddsTC <- DESeq(dds, full = mm.full, reduced = mm.reduced, test="LRT")
ddsTC <- DESeq(dds)
result <- results(ddsTC)
head(result)
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
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_100_DGE_nosalt24hs3_syl.pdf",width = 15,height = 15)

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
res_heat1h <- results(ddsTC, contrast=list("groupheat.1h", "groupCK.1h"))
head(res_heat1h)

#plot(res_heat1h$log2FoldChange,-log10(res_heat1h$padj))
#res_heat1h_lfc <- lfcShrink(ddsTC,coef = 20,res=res_heat1h)
#head(res_heat1h_lfc)
EnhancedVolcano(res_heat1h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Heat1H',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_heat1h_vocano_fc1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_heat2h <- results(ddsTC, contrast=list("groupheat.2h", "groupCK.2h"))
EnhancedVolcano(res_heat2h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Heat2H',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_heat2h_vocano_fc1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_heat5h <- results(ddsTC, contrast=list("groupheat.5h", "groupCK.5h"))              
EnhancedVolcano(res_heat5h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Heat5H',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_heat5h_vocano_fc1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_heat10h <- results(ddsTC, contrast=list("groupheat.10h", "groupCK.10h"))               
EnhancedVolcano(res_heat10h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Heat10H',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_heat10h_vocano_fc1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
res_heat24h <- results(ddsTC, contrast=list("groupheat.24h", "groupCK.24h"))
EnhancedVolcano(res_heat24h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Heat24H',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_heat24h_vocano_fc1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#Define the DEG as foldchanges >1 and adjustP-value <0.05
#heat1h
res_heat1h <- res_heat1h[!is.na(res_heat1h$padj),]
sigres_heat1h_up <- res_heat1h[(res_heat1h$padj<0.05 & res_heat1h$log2FoldChange>= 1),]
nb_heat1h_up <- nrow(sigres_heat1h_up)
#3160 #113 #358
write.table(x = sigres_heat1h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat1h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat1h_down <- res_heat1h[(res_heat1h$padj<0.05 & res_heat1h$log2FoldChange<= -1),]
nb_heat1h_down <- nrow(sigres_heat1h_down)
#3239 #0 #45 #2803 #actually is 3239
#xx<-res_heat1h[!is.na(res_heat1h$padj),]
#sigres_heat1h_down <-xx[(xx$log2FoldChange<=-1 & xx$padj <0.05),]
  
write.table(x = sigres_heat1h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat1h_down#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat1h_DEG <- res_heat1h[(res_heat1h$padj<0.05 & abs(res_heat1h$log2FoldChange) >= 1 ),]
nb_heat1h_DEG <- nrow(sigres_heat1h_DEG)
#6399
write.table(x = sigres_heat1h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat1h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

  
#heat2h
res_heat2h <- res_heat2h[!is.na(res_heat2h$padj),]
sigres_heat2h_up <- res_heat2h[(res_heat2h$padj< 0.05 & res_heat2h$log2FoldChange>= 1),]
nb_heat2h_up <- nrow(sigres_heat2h_up)
#2603 #37 #165
write.table(x = sigres_heat2h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat2h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat2h_down <- res_heat2h[(res_heat2h$padj<0.05 & res_heat2h$log2FoldChange<= -1),]
nb_heat2h_down <- nrow(sigres_heat2h_down)
#3077 #1 #62
write.table(x = sigres_heat2h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat2h_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat2h_DEG <- res_heat1h[(res_heat2h$padj<0.05 & abs(res_heat2h$log2FoldChange) >= 1 ),]
nb_heat2h_DEG <- nrow(sigres_heat2h_DEG)
#5680
write.table(x = sigres_heat2h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat2h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#heat5h
res_heat5h <- res_heat5h[!is.na(res_heat5h$padj),]
sigres_heat5h_up <- res_heat5h[(res_heat5h$padj<0.05 & res_heat5h$log2FoldChange>= 1),]
nb_heat5h_up <- nrow(sigres_heat5h_up)
#1865 #19
write.table(x = sigres_heat5h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat5h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_heat5h_down <- res_heat5h[(res_heat5h$padj<0.05 & res_heat5h$log2FoldChange<= -1),]
nb_heat5h_down <- nrow(sigres_heat5h_down)
#1784 #0
write.table(x = sigres_heat5h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat5h_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_heat5h_DEG <- res_heat5h[(res_heat5h$padj<0.05 & abs(res_heat5h$log2FoldChange) >= 1 ),]
nb_heat5h_DEG <- nrow(sigres_heat5h_DEG)
#3649
write.table(x = sigres_heat5h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat5h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


#heat10h
res_heat10h <- res_heat10h[!is.na(res_heat10h$padj),]
sigres_heat10h_up <- res_heat10h[(res_heat10h$padj < 0.05 & res_heat10h$log2FoldChange >= 1),]
nb_heat10h_up <- nrow(sigres_heat10h_up)
#1932 #10
write.table(x = sigres_heat10h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat10h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat10h_down <- res_heat10h[(res_heat10h$padj<0.05 & res_heat10h$log2FoldChange<= -1),]
nb_heat10h_down <- nrow(sigres_heat10h_down)
#1916 #0
write.table(x = sigres_heat10h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat10h_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat10h_DEG <- res_heat10h[(res_heat10h$padj<0.05 & abs(res_heat10h$log2FoldChange) >= 1 ),]
nb_heat10h_DEG <- nrow(sigres_heat10h_DEG)
#4225
write.table(x = sigres_heat10h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat10h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#heat24h
res_heat24h <- res_heat24h[!is.na(res_heat24h$padj),]
sigres_heat24h_up <- res_heat24h[(res_heat24h$padj<0.05 & res_heat24h$log2FoldChange>=1),]
nb_heat24h_up <- nrow(sigres_heat24h_up)
#3316 #9
write.table(x = sigres_heat24h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat24h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_heat24h_down <- res_heat24h[(res_heat24h$padj<0.05 & res_heat24h$log2FoldChange<= -1),]
nb_heat24h_down <- nrow(sigres_heat24h_down)
#2477 #
write.table(x = sigres_heat24h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_heat24h_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_heat24h_DEG <- res_heat24h[(res_heat24h$padj<0.05 & abs(res_heat24h$log2FoldChange) >= 1 ),]
nb_heat24h_DEG <- nrow(sigres_heat24h_DEG)
#5793
write.table(x = sigres_heat24h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat24h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

####Drought
res_drought1h <- results(ddsTC, contrast=list("groupdrought.1h", "groupCK.1h")) 
EnhancedVolcano(res_drought1h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Drought1H',
                xlim = c(-25, 25),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_drought1h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_drought2h <- results(ddsTC, contrast=list("groupdrought.2h", "groupCK.2h")) 
EnhancedVolcano(res_drought2h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Drought2H',
                xlim = c(-25, 25),
                FCcutoff=5,
                pCutoff = 10e-50,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_drought2h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_drought5h <- results(ddsTC, contrast=list("groupdrought.5h", "groupCK.5h"))
EnhancedVolcano(res_drought5h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Drought5H',
                xlim = c(-25, 25),
                FCcutoff=5,
                pCutoff = 10e-50,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_drought5h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_drought10h <- results(ddsTC, contrast=list("groupdrought.10h", "groupCK.10h")) 
EnhancedVolcano(res_drought10h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Drought10H',
                xlim = c(-25, 25),
                FCcutoff=5,
                pCutoff = 10e-50,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_drought10h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_drought24h <- results(ddsTC, contrast=list("groupdrought.24h", "groupCK.24h"))
EnhancedVolcano(res_drought24h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Drought24H',
                xlim = c(-25, 25),
                FCcutoff=5,
                pCutoff = 10e-50,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_drought24h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#drought
#heat1h
res_drought1h <- res_drought1h[!is.na(res_drought1h$padj),]
sigres_drought1h_up <- res_drought1h[(res_drought1h$padj<0.05 & res_drought1h$log2FoldChange> 1),]
nb_drought1h_up <- nrow(sigres_drought1h_up)
#2395 #0
sigres_drought1h_down <- res_drought1h[(res_drought1h$padj<0.05 & res_drought1h$log2FoldChange< -1),]
nb_drought1h_down <- nrow(sigres_drought1h_down)
#1361 #0
sigres_drought1h_DEG <- res_drought1h[(res_drought1h$padj<0.05 & abs(res_drought1h$log2FoldChange) >= 1 ),]
nb_drought1h_DEG <- nrow(sigres_drought1h_DEG)
#4225
write.table(x = sigres_drought1h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought1h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#heat2h
res_drought2h <- res_drought2h[!is.na(res_drought2h$padj),]
sigres_drought2h_up <- res_drought2h[(res_drought2h$padj<0.05 & res_drought2h$log2FoldChange> 1),]
nb_drought2h_up <- nrow(sigres_drought2h_up)
#3453 #2
sigres_drought2h_down <- res_drought2h[(res_drought2h$padj<0.05 & res_drought2h$log2FoldChange< -1),]
nb_drought2h_down <- nrow(sigres_drought2h_down)
#2544 #2
sigres_drought2h_DEG <- res_drought2h[(res_drought2h$padj<0.05 & abs(res_drought2h$log2FoldChange) >= 1 ),]
nb_drought2h_DEG <- nrow(sigres_drought2h_DEG)
#4225
write.table(x = sigres_drought2h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought2h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#heat5h
res_drought5h <- res_drought5h[!is.na(res_drought5h$padj),]

sigres_drought5h_up <- res_drought5h[(res_drought5h$padj<0.05 & res_drought5h$log2FoldChange> 1),]
nb_drought5h_up <- nrow(sigres_drought5h_up)
#5105 #38
sigres_drought5h_down <- res_drought5h[(res_drought5h$padj<0.05 & res_drought5h$log2FoldChange< -1),]
nb_drought5h_down <- nrow(sigres_drought5h_down)
#5017 #2
sigres_drought5h_DEG <- res_drought5h[(res_drought5h$padj<0.05 & abs(res_drought5h$log2FoldChange) >= 1 ),]
nb_drought5h_DEG <- nrow(sigres_drought5h_DEG)
#4225
write.table(x = sigres_drought5h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought5h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
#heat10h
res_drought10h <- res_drought10h[!is.na(res_drought10h$padj),]
sigres_drought10h_up <- res_drought10h[(res_drought10h$padj<0.05 & res_drought10h$log2FoldChange>1),]
nb_drought10h_up <- nrow(sigres_drought10h_up)
#5369 #31
sigres_drought10h_down <- res_drought10h[(res_drought10h$padj<0.05 & res_drought10h$log2FoldChange< -1),]
nb_drought10h_down <- nrow(sigres_drought10h_down)
#5014 #0
sigres_drought10h_DEG <- res_drought10h[(res_drought10h$padj<0.05 & abs(res_drought10h$log2FoldChange) >= 1 ),]
nb_drought10h_DEG <- nrow(sigres_drought10h_DEG)
#10383
write.table(x = sigres_drought10h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought10h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
            )

#heat24h
res_drought24h <- res_drought24h[!is.na(res_drought24h$padj),]

sigres_drought24h_up <- res_drought24h[(res_drought24h$padj<0.05 & res_drought24h$log2FoldChange>1),]
nb_drought24h_up <- nrow(sigres_drought24h_up)
#5853 #71
sigres_drought24h_down <- res_drought24h[(res_drought24h$padj<0.05 & res_drought24h$log2FoldChange< -1),]
nb_drought24h_down <- nrow(sigres_drought24h_down)
#5701 #4
sigres_drought24h_DEG <- res_drought24h[(res_drought24h$padj<0.05 & abs(res_drought24h$log2FoldChange) >= 1 ),]
nb_drought24h_DEG <- nrow(sigres_drought24h_DEG)
#11554

write.table(x = sigres_drought24h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought24h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

####salt
res_salt1h <- results(ddsTC, contrast=list("groupsalt.1h", "groupCK.1h"))
EnhancedVolcano(res_salt1h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Salt1H',
                xlim = c(-25, 25),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_salt1h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_salt2h <- results(ddsTC, contrast=list("groupsalt.2h", "groupCK.2h")) 
EnhancedVolcano(res_salt2h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Salt2H',
                xlim = c(-25, 25),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_salt2h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_salt5h <- results(ddsTC, contrast=list("groupsalt.5h", "groupCK.5h")) 
EnhancedVolcano(res_salt5h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Salt5H',
                xlim = c(-25, 25),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_salt5h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
res_salt10h <- results(ddsTC, contrast=list("groupsalt.10h", "groupCK.10h")) 
EnhancedVolcano(res_salt10h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Salt10H',
                xlim = c(-25, 25),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_salt10h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_salt24h <- results(ddsTC, contrast=list("groupsalt.24h", "groupCK.24h")) 
EnhancedVolcano(res_salt24h,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Salt24H',
                xlim = c(-25, 25),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/res_salt24h_vocano#.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

#na.omit(res_heat1h)
#heat1h
res_salt1h <- res_salt1h[!is.na(res_salt1h$padj),]
sigres_salt1h_up <- res_salt1h[(res_salt1h$padj<0.05 & res_salt1h$log2FoldChange>=1),]
nb_salt1h_up <- nrow(sigres_salt1h_up)
#1460 #0
sigres_salt1h_down <- res_salt1h[(res_salt1h$padj<0.05 & res_salt1h$log2FoldChange <= -1),]
nb_salt1h_down <- nrow(sigres_salt1h_down)
#188
sigres_salt1h_DEG <- res_salt1h[(res_salt1h$padj<0.05 & abs(res_salt1h$log2FoldChange) >= 1 ),]
nb_salt1h_DEG <- nrow(sigres_salt1h_DEG)
#1648

write.table(x = sigres_salt1h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt1h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#heat2h
res_salt2h <- res_salt2h[!is.na(res_salt2h$padj),]
sigres_salt2h_up <- res_salt2h[(res_salt2h$padj<0.05 & res_salt2h$log2FoldChange>1),]
nb_salt2h_up <- nrow(sigres_salt2h_up)
#594
sigres_salt2h_down <- res_salt2h[(res_salt2h$padj<0.05 & res_salt2h$log2FoldChange< -1),]
nb_salt2h_down <- nrow(sigres_salt2h_down)
#136
sigres_salt2h_DEG <- res_salt1h[(res_salt2h$padj<0.05 & abs(res_salt2h$log2FoldChange) >= 1 ),]
nb_salt2h_DEG <- nrow(sigres_salt2h_DEG)
#1648

write.table(x = sigres_salt2h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt2h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


#heat5h
res_salt5h <- res_salt5h[!is.na(res_salt5h$padj),]
sigres_salt5h_up <- res_salt5h[(res_salt5h$padj<0.05 & res_salt5h$log2FoldChange> 1),]
nb_salt5h_up <- nrow(sigres_salt5h_up)
#301
sigres_salt5h_down <- res_salt5h[(res_salt5h$padj<0.05 & res_salt5h$log2FoldChange< -1),]
nb_salt5h_down <- nrow(sigres_salt5h_down)
#71
sigres_salt5h_DEG <- res_salt5h[(res_salt5h$padj<0.05 & abs(res_salt5h$log2FoldChange) >= 1 ),]
nb_salt5h_DEG <- nrow(sigres_salt5h_DEG)
#1648

write.table(x = sigres_salt5h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt5h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#heat10h
res_salt10h <- res_salt10h[!is.na(res_salt10h$padj),]

sigres_salt10h_up <- res_salt10h[(res_salt10h$padj<0.05 & res_salt10h$log2FoldChange>1),]
nb_salt10h_up <- nrow(sigres_salt10h_up)
#491
sigres_salt10h_down <- res_salt10h[(res_salt10h$padj<0.05 & res_salt10h$log2FoldChange< -1),]
nb_salt10h_down <- nrow(sigres_salt10h_down)
#292
sigres_salt10h_DEG <- res_salt10h[(res_salt10h$padj<0.05 & abs(res_salt10h$log2FoldChange) >= 1 ),]
nb_salt10h_DEG <- nrow(sigres_salt10h_DEG)
#783

write.table(x = sigres_salt10h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt10h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#heat24h
res_salt24h <- res_salt24h[!is.na(res_salt24h$padj),]

sigres_salt24h_up <- res_salt24h[(res_salt24h$padj<0.05 & res_salt24h$log2FoldChange>1),]
nb_salt24h_up <- nrow(sigres_salt24h_up)
#793
sigres_salt24h_down <- res_salt24h[(res_salt24h$padj<0.05 & res_salt24h$log2FoldChange< -1),]
nb_salt24h_down <- nrow(sigres_salt24h_down)
#633
sigres_salt24h_DEG <- res_salt24h[(res_salt24h$padj<0.05 & abs(res_salt24h$log2FoldChange) >= 1 ),]
nb_salt24h_DEG <- nrow(sigres_salt24h_DEG)
#1426

write.table(x = sigres_salt24h_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt24h_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)



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

pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat_DEG_nb_nosalt24hs3_syl.pdf",width = 8,height = 8)
p <- ggplot(data=heat, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,6000)
  
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

pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought_DEG_nb_nosalt_24hs3_syl.pdf",width = 8,height = 8)
p <- ggplot(data=drought, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,6000)

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

pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt_DEG_nb_no_salt_24h_s3_syl.pdf",width = 8,height = 8)
p <- ggplot(data=salt, aes(x=Time, y=Number, fill=Categories)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,6000)

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
#5987
#downregulating:
heat1_2_down <- union(rownames(sigres_heat1h_down),rownames(sigres_heat2h_down))
heat1_5_down <- union(heat1_2_down,rownames(sigres_heat5h_down))
heat1_10_down <- union(heat1_5_down,rownames(sigres_heat10h_down))
heat1_24_down <- union(heat1_10_down,rownames(sigres_heat24h_down))
length(heat1_24_down)
#5603
####Drought:
#Take the union of the upregulating genes under drought tretment:
#upregulating:
drought1_2_up <- union(rownames(sigres_drought1h_up),rownames(sigres_drought2h_up))
drought1_5_up <- union(drought1_2_up,rownames(sigres_drought5h_up))
drought1_10_up <- union(drought1_5_up,rownames(sigres_drought10h_up))
drought1_24_up <- union(drought1_10_up,rownames(sigres_drought24h_up))
length(drought1_24_up)
#7384
#downregulating:
drought1_2_down <- union(rownames(sigres_drought1h_down),rownames(sigres_drought2h_down))
drought1_5_down <- union(drought1_2_down,rownames(sigres_drought5h_down))
drought1_10_down <- union(drought1_5_down,rownames(sigres_drought10h_down))
drought1_24_down <- union(drought1_10_down,rownames(sigres_drought24h_down))
length(drought1_24_down)
#7671

#Salt:
#Take the union of the upregulating genes under salt tretment:
#upregulating:
salt1_2_up <- union(rownames(sigres_salt1h_up),rownames(sigres_salt2h_up))
salt1_5_up <- union(salt1_2_up,rownames(sigres_salt5h_up))
salt1_10_up <- union(salt1_5_up,rownames(sigres_salt10h_up))
salt1_24_up <- union(salt1_10_up,rownames(sigres_salt24h_up))
length(salt1_24_up)
#2547
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
setwd("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/")
#dev.off()
venn.diagram(
  x = list(heat1_24_up, drought1_24_up, salt1_24_up),
  category.names = c("Heat" , "Drought " , "Salt"),
  filename = 'up-DGE-abiotic_syl.png',
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
  filename = 'down-DGE-abiotic_syl.png',
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
#Up-regulated
#2785+1925+3286+1041+1132+145+1041+229=11584
#down-regulated
#2220+2647+4250+580+156+194+78=10125
#heat
#2785+1925+1132+145=5987
#drought
#1925+3286+1132+1041=7384
#salt
#1132+145+1041+229=2547

#down
#heat
#2220+2647+580+156=5603
#drought
#2647+4250+580+194=7671
#sailt
#156+580+194+78=1008

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
#Write the data into the excel

write.table(
  combined,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/DGE_shoot_abiotic_sylvaticum.txt",
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
             
