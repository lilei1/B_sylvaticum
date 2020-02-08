#11-26-2019 Berkeley by Li Lei
#This is for normalization for the RNAseq data subset of the abioticstress according to Avinash's sugestions.
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("rlang")
library(edgeR)
#packageVersion("rlang")
counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
head(counts)

head <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tmp_counts_head.txt",header = F)
colnames(counts) <- head$V2 #change the headers
abiotic <- subset(counts,select = c(1:60))#take the subset of the data
head(abiotic)

d0 <- DGEList(abiotic)
d0 <- calcNormFactors(d0)
d0

cutoff <- 13
#test_cpm <- cpm(d0,normalized.lib.sizes = TRUE,prior.count = 0.25,log = TRUE)
#head(test_cpm)
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d)
head(d)
snames <- colnames(abiotic)

treatments <- sapply(strsplit(snames, "_",fixed = TRUE),function(x) (x[2]))
time <- sapply(strsplit(snames, "_",fixed = TRUE),function(x) (x[3]))
time <- sapply(strsplit(time, "h",fixed = TRUE),function(x) (x[1]))#use h to split the name

group <- interaction(treatments, time)
group
####
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_repli.pdf",width = 25,height = 25)
plotMDS(d, col = as.numeric(group))
dev.off()


mm <- model.matrix(~0 + group)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/abiotic_voom.pdf",width = 25,height = 25)
y <- voom(d, mm, plot = T)
dev.off()
voom_count <- data.frame(y$E)
#y <- voom(d0, mm, plot = T)
head(voom_count)
write.table(voom_count, file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/voom_normalized_count.txt", sep="\t",quote = F) 
#fit <- lmFit(y, mm)
#fit <- eBayes(fit)
#head(coef(fit))
#p.adjusted <- p.adjust(fit$p.value[,2], method="fdr") 
#results_limma <- cbind(fit$coeff, fit$p.value[,2], p.adjusted) 
#colnames(results_limma) <- c("av_expr", "2LogFC", "pvalue", "adjusted_pvalue") 
#head(results_limma)
