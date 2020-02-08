source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("rlang")
biocLite("DESeq2")
library(edgeR)
library(DESeq2)
#packageVersion("rlang")
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)

head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
abiotic <- subset(counts,select = c(1:60))#take the subset of the data
head(abiotic)
#This is for normalization for the RNAseq data subset of the abioticstress according to Avinash's sugestions.
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("rlang")
biocLite("DESeq2")
library(edgeR)
#packageVersion("rlang")
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#36927

head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
counts <- subset(counts,select = c(1:60))#take the subset of the data

meta.data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/meta_abio.txt",header = T)
head(meta.data)

head(counts)
count.cutoff = 5 # 10
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]

nrow(counts)
#25072
#CPM#
count.cutoff = 1
bioreplicates.cutoff = 3
normalized.counts <- cpm(counts,lib.size = NULL)
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#21821
counts <- counts[keep, ]

nrow(counts)
##25072
row.names(normalized.counts)
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#[1] 21821
counts[21822,]
grep("Brasy3G282900", rownames(counts))
#21822/36927
## VST instead of voom
#counts[is.na(counts)] <- 0
dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = as.data.frame(meta.data), design = ~group )
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#row.names(expr.vst)
#nrow(expr.vst)
#21821
#grep("Brasy3G282900", rownames(expr.vst))
write.table(expr.vst, file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_normalized_count_#1.txt",sep="\t",quote = F) 
#fit <- lmFit(y, mm)
#fit <- eBayes(fit)
#head(coef(fit))
#p.adjusted <- p.adjust(fit$p.value[,2], method="fdr") 
#results_limma <- cbind(fit$coeff, fit$p.value[,2], p.adjusted) 
#colnames(results_limma) <- c("av_expr", "2LogFC", "pvalue", "adjusted_pvalue") 
#head(results_limma)
