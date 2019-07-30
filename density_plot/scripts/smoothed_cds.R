#!/usr/bin/env Rscript
#   Script to calculate Lowess-smoothed estimates of CDS density in B.Sylvaticum.
# Written by Li Lei, July 30th, 2019 in Walnut Creek, CA.

#Define a function to read file
readFile <- function(filename) {
  cdsData <- read.table(file = filename,
                           header = TRUE) # include column names
  return(cdsData)
}

#   define a function to return the smoothed values
smooth <- function(chrom, mapdata, f, delta){
  map_subset <- mapdata[mapdata$Chromosome == chrom,]
  midpoints <- (map_subset$LeftBP + map_subset$RightBP)/2
  smoothed <- lowess(map_subset$CDS.Mb~midpoints, f=f, delta=delta)
  smoothed$y[smoothed$y < 0] <- 0
  return(smoothed)
}

#setwd("/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts")
#cmmb <- read.table("chr1_Bsylvaticum_490_v1.1_cds_rpt1.txt", header=T)

#head(cmmb)

#define new interval:
interval <- function(cmmb){
chromosomes <- unique(cmmb$Chromosome)
smoothed_map <- sapply(chromosomes, smooth, cmmb, f=0.02, delta=3000000)

#   Associate the smoothed values with the original map intervals
new_intervals <- sapply(
  seq_along(chromosomes),
  function(chrom) {
    map_subset <- cmmb[cmmb$Chromosome==chromosomes[chrom], ]
    map_subset$CDS.Mb <- smoothed_map["y", chrom]$y
    return(map_subset)
  })

df_new_intervals <- data.frame(
  Chromosome=as.character(unlist(new_intervals["Chromosome",])),
  LeftBP=as.numeric(unlist(new_intervals["LeftBP",])),
  RightBP=as.numeric(unlist(new_intervals["RightBP",])),
  Smoothed_cMMb=as.numeric(unlist(new_intervals["CDS.Mb",]))
)
return(df_new_intervals)
}
#   Save the new file!
writeFile <- function (df_new_intervals,output){
  write.table(
  df_new_intervals,
  file=output,
  quote=F,
  row.names=F,
  sep="\t")
}

#   Driver function
main <- function() {
  #   Take command line arguments
  #   Stores arguments into a vector
  args <- commandArgs(trailingOnly = TRUE)
  #   User provided arguments
  genoData <- args[1] # 1) genotype data frame that has SNPs sorted
  #inDir <- args[2] # 2) Name of whole chr used for heatmap
  #inName <- args[2]
  outDir <- args[2] # 5) where do our output files go?
  outName <- args[3] # 7) Do we want to include SNP names in our plot? (include/exclude)
  
  #infile = paste0(inDir, "/", inName)
  outfile = paste0(outDir, "/", outName, ".txt")
  #   Read in genotype data and SNP_BAC.txt file
  genoFile <- readFile(filename=genoData)
  new_intervals <- interval(cmmb=genoFile)
  #   Apply function over every SNP in inverted region
  finalfile <- writeFile(df_new_intervals=new_intervals,output=outfile)
}

#   Run program
main()
