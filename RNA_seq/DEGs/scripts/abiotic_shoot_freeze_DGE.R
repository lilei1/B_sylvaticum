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

#head <- read.delim("/global/projectb/scratch/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
#colnames(counts) <- head$V2 #change the headers
#Just focus on the shoot along the time point
counts <- subset(counts,select = c(1:60))#take the subset of the data
#counts$Shoot_salt_24h.s3 <- NULL
#meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
meta.data <- read.delim("/global/projectb/scratch/llei2019/Brachypodium/Sylvaticum/RNAseq/coldata_abiotic_shoot.txt",header = T)
head(meta.data)
nrow(meta.data)
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
count.cutoff = 3 # 60% of the data were above 3
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#[1] 26421
normalized.counts <- cpm(counts)
head(normalized.counts)
normalized.counts.re <- data.frame(normalized.counts)
head(normalized.counts.re)
normalized.counts_long <- gather(normalized.counts.re,code,counts)
quantile(normalized.counts_long$counts,probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#7%           8%           9%          10%          11%          12%          13% 
#0.000000e+00 0.000000e+00 5.168025e-02 5.837050e-02 6.294471e-02 7.050105e-02 9.484186e-02 
#14%          15%          16%          17%          18%          19%          20% 
#1.134373e-01 1.319454e-01 1.549899e-01 1.751115e-01 2.015789e-01 2.308456e-01 2.584012e-01 
#21%          22%          23%          24%          25%          26%          27% 
#2.933894e-01 3.345268e-01 3.723011e-01 4.204855e-01 4.742093e-01 5.253345e-01 5.837502e-01 
#28%          29%          30%          31%          32%          33%          34% 
#6.502081e-01 7.242551e-01 8.059097e-01 8.863553e-01 9.853936e-01 1.084930e+00 1.198233e+00 
#35%          36%          37%          38%          39%          40%          41% 
#1.314175e+00 1.449616e+00 1.588921e+00 1.747017e+00 1.904435e+00 2.084500e+00 2.273768e+00 
#42%          43%          44%          45%          46%          47%          48% 
#2.475294e+00 2.703537e+00 2.935814e+00 3.182274e+00 3.452616e+00 3.735666e+00 4.034104e+00 
#49%          50%          51%          52%          53%          54%          55% 
#4.361200e+00 4.703509e+00 5.068935e+00 5.442631e+00 5.850763e+00 6.282443e+00 6.736300e+00 
#56%          57%          58%          59%          60%          61%          62% 
#7.208094e+00 7.715582e+00 8.246659e+00 8.809321e+00 9.397890e+00 1.002064e+01 1.068758e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.138758e+01 1.212329e+01 1.290105e+01 1.371986e+01 1.459372e+01 1.551036e+01 1.648214e+01 
#70%          71%          72%          73%          74%          75%          76% 
#1.750283e+01 1.858619e+01 1.974997e+01 2.097001e+01 2.229044e+01 2.369034e+01 2.520973e+01 
#77%          78%          79%          80%          81%          82%          83% 
#2.681399e+01 2.852955e+01 3.039431e+01 3.237739e+01 3.454449e+01 3.685717e+01 3.934653e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.212751e+01 4.520064e+01 4.862248e+01 5.244484e+01 5.678486e+01 6.179615e+01 6.769817e+01 
#91%          92%          93%          94%          95%          96%          97% 
#7.457202e+01 8.312925e+01 9.359878e+01 1.068847e+02 1.245855e+02 1.504111e+02 1.912071e+02 
#98%          99%         100% 
#2.634619e+02 4.568648e+02 2.980066e+04
head(normalized.counts_long)
count.cutoff = 1 #86% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#[1] 21821
nrow(counts)
#26421
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#21821
grep("Brasy1G039900", rownames(counts))#check if this CBF got filtered, CBF alos for drought
#[1] 255
ncol(counts)
#60
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
#dev.off()
pdf(file = "/global/projectb/scratch/llei2019/Brachypodium/Sylvaticum/RNAseq/heatmap_100_DGE_vst_syl.pdf",width = 15,height = 15)

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

head <- read.delim("/global/projectb/scratch/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
counts <- subset(counts,select = c(1:60))#take the subset of the data
counts$Shoot_salt_24h.s3 <- NULL
#meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
meta.data <- read.delim("/global/projectb/scratch/llei2019/Brachypodium/Sylvaticum/RNAseq/coldata_abiotic_shoot_nosalt24hs3.txt",header = T)
head(meta.data)
group <- paste0(meta.data$treatment, ".", meta.data$time)
####do filtering with raw count
count.cutoff = 3 # 10 #above 60% of the datapoint is greater than 3;
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#26386
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
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#7%           8%           9%          10%          11%          12%          13% 
#0.000000e+00 0.000000e+00 5.235372e-02 5.910984e-02 6.344763e-02 7.334741e-02 9.643267e-02 
#14%          15%          16%          17%          18%          19%          20% 
#1.181810e-01 1.363530e-01 1.573112e-01 1.783206e-01 2.050308e-01 2.371060e-01 2.637600e-01 
#21%          22%          23%          24%          25%          26%          27% 
#2.997233e-01 3.391749e-01 3.793678e-01 4.245458e-01 4.821634e-01 5.347530e-01 5.944021e-01 
#28%          29%          30%          31%          32%          33%          34% 
#6.638936e-01 7.364594e-01 8.182815e-01 9.036288e-01 9.963042e-01 1.102318e+00 1.209632e+00 
#35%          36%          37%          38%          39%          40%          41% 
#1.332400e+00 1.468238e+00 1.612313e+00 1.763523e+00 1.929351e+00 2.110080e+00 2.301746e+00 
#42%          43%          44%          45%          46%          47%          48% 
#2.507249e+00 2.727061e+00 2.967300e+00 3.223724e+00 3.491775e+00 3.776889e+00 4.087116e+00 
#49%          50%          51%          52%          53%          54%          55% 
#4.405009e+00 4.748449e+00 5.111858e+00 5.498897e+00 5.902391e+00 6.335618e+00 6.793728e+00 
#56%          57%          58%          59%          60%          61%          62% 
#7.261393e+00 7.778986e+00 8.307065e+00 8.863575e+00 9.460074e+00 1.008167e+01 1.075171e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.145211e+01 1.219117e+01 1.296912e+01 1.378981e+01 1.466035e+01 1.557704e+01 1.654378e+01 
#70%          71%          72%          73%          74%          75%          76% 
#1.756642e+01 1.865082e+01 1.982197e+01 2.103486e+01 2.234953e+01 2.375831e+01 2.526756e+01 
#77%          78%          79%          80%          81%          82%          83% 
#2.686494e+01 2.858482e+01 3.044196e+01 3.242007e+01 3.458186e+01 3.688242e+01 3.936675e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.213333e+01 4.520423e+01 4.862209e+01 5.243141e+01 5.675808e+01 6.175047e+01 6.762812e+01 
#91%          92%          93%          94%          95%          96%          97% 
#7.448751e+01 8.300526e+01 9.340574e+01 1.066931e+02 1.243500e+02 1.500766e+02 1.908709e+02 
#98%          99%         100% 
#2.629703e+02 4.562462e+02 2.980068e+04  
count.cutoff = 1 #86% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#21794
nrow(counts)
#26386
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#21794
#counts[21822,]
grep("Brasy3G282900", rownames(counts))#check if this gene hsp gene got filtered
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


nrow(dds)
#21794
#nrow(dds2)
#21790

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
write.table(x = sigres_drought1h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought1h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_drought1h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought1h_down",
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
write.table(x = sigres_drought2h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought2h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_drought2h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought2h_down.txt",
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
write.table(x = sigres_drought5h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought5h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_drought5h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought5h_down.txt",
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
write.table(x = sigres_drought10h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought10h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_drought10h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought10h_down.txt",
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
write.table(x = sigres_drought24h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought24h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_drought24h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_drought24h_down.txt",
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
write.table(x = sigres_salt1h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt1h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_salt1h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt1h_down.txt",
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
write.table(x = sigres_salt2h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt2h_uptxt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_salt2h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt2h_down.txt",
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
write.table(x = sigres_salt5h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt5h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_salt5h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt5h_down.txt",
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
write.table(x = sigres_salt10h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt10h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_salt10h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt10h_down.txt",
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
write.table(x = sigres_salt24h_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt24h_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = sigres_salt24h_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_salt24h_down.txt",
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

pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat_DEG_nb_nosalt24hs3_syl_combined.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought_DEG_nb_nosalt_24hs3_syl_combined.pdf",width = 8,height = 8)
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

pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt_DEG_nb_no_salt_24h_s3_syl_combined.pdf",width = 8,height = 8)
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
# 2547
#downregulating:
salt1_2_down <- union(rownames(sigres_salt1h_down),rownames(sigres_salt2h_down))
salt1_5_down <- union(salt1_2_down,rownames(sigres_salt5h_down))
salt1_10_down <- union(salt1_5_down,rownames(sigres_salt10h_down))
salt1_24_down <- union(salt1_10_down,rownames(sigres_salt24h_down))
length(salt1_24_down)
#1008

write.table(x = salt1_24_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt1_24_up#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = drought1_24_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought1_24_up#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = heat1_24_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat1_24_up#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

write.table(x = heat1_24_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heat1_24_down#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = drought1_24_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/drought1_24_down#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
write.table(x = salt1_24_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/salt1_24_down#.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
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
  x = list(drought1_24_up,heat1_24_up,  salt1_24_up),
  category.names = c( "Drought " , "Heat" ,"Salt"),
  filename = 'up-DGE-abiotic_syl#.png',
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
  x = list( drought1_24_down, heat1_24_down,salt1_24_down),
  category.names = c("Drought " , "Heat" , "Salt"),
  filename = 'down-DGE-abiotic_syl#.png',
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

salt_DEG_gene_list <- union(union(union(union(rownames(sigres_salt1h_DEG),rownames(sigres_salt2h_DEG)),rownames(sigres_salt5h_DEG)),rownames(sigres_salt10h_DEG)),rownames(sigres_salt24h_DEG))
length(salt_DEG_gene_list)
#3512
drought_DEG_gene_list <- union(union(union(union(rownames(sigres_drought1h_DEG),rownames(sigres_drought2h_DEG)),rownames(sigres_drought5h_DEG)),rownames(sigres_drought10h_DEG)),rownames(sigres_drought24h_DEG))
length(drought_DEG_gene_list)
#14694

heat_DEG_gene_list <- union(union(union(union(rownames(sigres_heat1h_DEG),rownames(sigres_heat2h_DEG)),rownames(sigres_heat5h_DEG)),rownames(sigres_heat10h_DEG)),rownames(sigres_heat24h_DEG))
length(heat_DEG_gene_list)
#11217
top_gene_list <- union(union(salt_DEG_gene_list,drought_DEG_gene_list),heat_DEG_gene_list)
length(top_gene_list)
#17184
write.table(x = top_gene_list,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/all_sig_DEGs_across_stresses_time_gene.list.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#plot the genes with logfolderchanges >1
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_17184_DGE_shoot_abiotic.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
head(assay(rld))
palette <- colorRampPalette(c("red","white","blue"))(256)
head(assay(rld))
reorder_rld <- assay(rld)[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
head(reorder_rld)
heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='none', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", drought="#ffff00",heat="#ff8000", salt="#0066ff")[
             colData(rld)$treatment ] )
dev.off()

#up-regulating
3286+1925+2785+1132+1041+145+229
#10543

#Down regulating
4250+2647+2220+580+194+156+78
#10125

#drought-up
3286+1925+1132+1041
#7384

#Heat-up
1925+2785+1132+145
#6987
#salt
1041+1132+145+229
#2547

##drought-down
4250+2647+580+194
#7671

#heat_down
2647+2220+580
#5447

#salt -down
194+580+56+78
#908
