file <- read.delim(file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/scripts/chr1_Bsylvaticum_490_v1.1_cds_rpt1.txt",header=F,sep = "\t")
head(file)
plot(file$V4/1000000,file$V5/1000,xlab = "Ch1(Mb)",ylab="CDS length(Kb)",pch=20)
y <- file$V5/1000
x <- file$V4/1000000
lo <- loess(y ~ x)
lines(predict(lo), col='red', lwd=2)
library("ggplot2")
qplot(x,y, geom='smooth', span =0.04)
