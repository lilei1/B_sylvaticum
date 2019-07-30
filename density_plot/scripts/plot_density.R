#!/usr/bin/env Rscript
#   Script to write a fancy plot for the smoothied CDS density file
# Written by Li Lei, July 30th, 2019 in Walnut Creek, CA.

#load GGplot2
#install.packages("ggplot2")
library(ggplot2)

#Define a function to read file
readFile <- function(filename) {
  cdsData <- read.table(file = filename,
                        header = TRUE) # include column names
  return(cdsData)
}


#rec_rate <- read.table("/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts/chr1_Bsylvaticum_490_v1.1_cds_rpt1.smooth.txt", header=T)
#head(rec_rate)

#Define a plot function:
plot <- function(rec_rate, outDir,outName){
pdf(file = paste0(outDir, "/", outName, ".pdf"),width=16.00, height=8.67)
ggplot(rec_rate) +
  geom_line(aes(x=(LeftBP+RightBP)/2000000, y=Smoothed_cMMb/1000), data=rec_rate, color="#9900ff", size=1, alpha=1) +
  facet_grid(Chromosome~.) +
  scale_y_continuous(limits=c(0, 100), breaks=seq(0,100,by=5)) +
  scale_x_continuous(limits=c(0, 50), breaks=seq(0, 50, by=5)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    strip.text.y=element_text(size=10, colour="black", angle=0),
    axis.ticks.y=element_blank()) +
  labs(y="", x="Physical Position (Mb)")

dev.off()
}

#   Driver function
main <- function() {
  #   Take command line arguments
  #   Stores arguments into a vector
  args <- commandArgs(trailingOnly = TRUE)
  #   User provided arguments
  genoData <- args[1] # 1) genotype data frame that has SNPs sorted
  outDir <- args[2] # 5) where do our output files go?
  outName <- args[3] # 7) Do we want to include SNP names in our plot? (include/exclude)
  input <- readFile(filename =genoData )
  #genoData <- "/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts/chr1_Bsylvaticum_490_v1.1_cds_rpt1.smooth.txt"
  #outDir ="/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts"
  #outName ="test"
  cds_plot <- plot(rec_rate = input, outDir = outDir, outName = outName)

}

#   Run program
main()
