#!/usr/bin/env Rscript
#   Script to do normalization for RNAseq data with DEseq.
# Written by Li Lei, 02-05, 2020 in Berkley, CA.
#getwd()
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("rlang")
#biocLite("DESeq2")
BiocManager::install(c("edgeR","rlang","DESeq2"))
library(edgeR)
library(DESeq2)

#Define a function to read raw count file
readCountFile <- function(filename) {
  rawcounts <- read.delim(file = filename,
                        row.names = 1,
                        header = T) # include column names
  return(rawcounts)
}

#Define a function to read lsmetadata file
readMetaFile <- function(filename) {
  meta_data <- read.delim(file = filename,
                          header = T) # include column names
  return(meta_data)
}

#define a function for writing file
writeOutFile <- function(normalCount, outputName) {
  write.table(x = normalCount,
              file = outputName,
              quote = FALSE,
              sep = "\t",
              eol = "\n",
              col.names = TRUE,
              row.names = FALSE)
}



#counts <- read.delim("/global/projectb/scratch/vrsingan/SeanGordon_Sylvaticum/gene_counts/counts.txt", row.names = 1)
#head(counts)
#nrow(counts)
#36927
#meta_data <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts/Meta_data_syl.txt",header = T)
#head(meta_data)
#nrow(abiotic_shoot)
#72
#abiotic_shoot <- meta_data[(meta_data$tissue == "shoot" & (meta_data$experiment == "timecourse" | meta_data$experiment == "abiotic")),]
#head(abiotic_shoot)

#Normalization function!!
norm_exp <- function(metaData,counts){
  #extract the column I needed!!!
  col.num <- which(metaData$code %in% colnames(counts))
  counts <- counts[, col.num]
  mata <- metaData[col.num,]
  
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
  #nrow(normalized.counts)
  #24110
  counts <- counts[keep, ]
  
  #nrow(counts)
  ##24110
  #row.names(normalized.counts)
  counts <- counts[row.names(normalized.counts),]
  #nrow(counts)
  #[1] 24110
  #counts[21822,]
  #grep("Brasy3G282900", rownames(counts))#check if this gene got filtered
  #21822/36927
  ## VST instead of voom
  #counts[is.na(counts)] <- 0
  #creat a group 
  mata$group <- paste(mata$treatment, mata$timepoint, sep='_')
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = as.data.frame(mata), design = ~group )
  expr.vst <- assay(DESeq2::vst(dds))
  #head(expr.vst)
  #str(expr.vst)
  expr.vst <- as.data.frame(expr.vst)
  return(expr.vst)
}

#row.names(expr.vst)
#nrow(expr.vst)
#21821
#grep("Brasy3G282900", rownames(expr.vst))
#write.table(expr.vst, file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/VST_normalized_count_abiotic.txt",sep="\t",quote = F) 

#   Function to run the script
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  rawCounts <- args[1]
  meta_data <- args[2]
  output <- args[3]
  counts <- readCountFile(rawCounts)
  metaData <- readMetaFile(meta_data)
  vst_norm_exp <- norm_exp(metaData,counts)
  writeOutFile(vst_norm_exp, output)
}

main() # Run program

