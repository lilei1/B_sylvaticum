#!/usr/bin/env Rscript
#   Script to do normalization for RNAseq data with DEseq.
# Written by Li Lei, 01-24, 2020 in Berkley, CA.
library(edgeR)
library(DESeq2)
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)
nrow(counts)
#36927
meta_data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts/Meta_data_syl.txt",header = T)
head(meta_data)
nrow(abiotic_shoot)
#72
abiotic_shoot <- meta_data[(meta_data$tissue == "shoot" & (meta_data$experiment == "timecourse" | meta_data$experiment == "abiotic")),]
head(abiotic_shoot)
#extract the column I needed!!!
col.num <- which(abiotic_shoot$code %in% colnames(counts))
counts <- counts[, col.num]

ncol(counts)
head(counts)

count.cutoff = 5 # 10
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]

nrow(counts)
#27323
#CPM#
count.cutoff = 1
bioreplicates.cutoff = 3
normalized.counts <- cpm(counts)
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#24110
counts <- counts[keep, ]

nrow(counts)
##24110
row.names(normalized.counts)
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#[1] 24110
counts[21822,]
grep("Brasy3G282900", rownames(counts))#check if this gene got filtered
#21822/36927
## VST instead of voom
#counts[is.na(counts)] <- 0
#creat a group 
abiotic_shoot$group <- paste(abiotic_shoot$treatment, abiotic_shoot$timepoint, sep='_')

dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = as.data.frame(abiotic_shoot), design = ~group )
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#row.names(expr.vst)
#nrow(expr.vst)
#21821
#grep("Brasy3G282900", rownames(expr.vst))
write.table(expr.vst, file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_normalized_count_abiotic.txt",sep="\t",quote = F) 
#fit <- lmFit(y, mm)
#fit <- eBayes(fit)
#head(coef(fit))
#p.adjusted <- p.adjust(fit$p.value[,2], method="fdr") 
#results_limma <- cbind(fit$coeff, fit$p.value[,2], p.adjusted) 
#colnames(results_limma) <- c("av_expr", "2LogFC", "pvalue", "adjusted_pvalue") 
#head(results_limma)

