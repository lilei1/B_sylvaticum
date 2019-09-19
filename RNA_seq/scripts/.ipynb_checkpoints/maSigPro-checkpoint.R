if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("maSigPro")
library(maSigPro)

data(NBdata)
head(NBdata)
data(NBdesign)
head(NBdesign)

###################
data(data.abiotic)
data(edesign.abiotic)
edesign.abiotic
rownames(edesign.abiotic)
colnames(edesign.abiotic)
head(data.abiotic)
design <- make.design.matrix(edesign.abiotic, degree = 4)
head(design)
design$groups.vector
library(MASS)

fit <- p.vector(data.abiotic, design, counts=TRUE)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
names(sigs)
names(sigs$sig.genes)
names(sigs$sig.genes$ColdvsControl)
suma2Venn(sigs$summary[, c(2:4)])
suma2Venn(sigs$summary[, c(1:4)])
see.genes(sigs$sig.genes$ColdvsControl, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 10)
#install.packages("XQuartz")
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)
PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T, dis = design$dis, groups.vector = design$groups.vector)
?see.genes
data(NBdata)
##############
edesign.abiotic <- data.frame(read.delim(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/design.abiotic",sep="\t", header=T))
rownames(edesign.abiotic) <- edesign.abiotic$X
edesign.abiotic <- edesign.abiotic[,(2:7)]
edesign.abiotic 

data <- read.delim(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/counts.txt", sep="\t", header=T)
data.abiotic <- data.frame(data[,c(1:61)])
rownames(data.abiotic) <- data.abiotic$GeneID
data.abiotic <- data.abiotic[,c(2:61)]
head(data.abiotic)
data.abiotic$median <- apply(data.abiotic, 1, median) 
subset <- data.abiotic[(data.abiotic$median<=1340.70 & data.abiotic$median>=10),]
head(subset)
sub.abiotic <-subset[,c(1:60)]
head(sub.abiotic)
design <- make.design.matrix(edesign.abiotic, degree = 4)

#nrow(data.abiotic)
#make the first column as row name

edesign.abiotic

quantile(data.abiotic$median,probs = seq(0,1,0.01))
#0%        1%        2%        3%        4%        5%        6%        7%        8%        9%       10%       11%       12%       13%       14%       15%       16% 
#0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00 
#17%       18%       19%       20%       21%       22%       23%       24%       25%       26%       27%       28%       29%       30%       31%       32%       33% 
#0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.50 
#34%       35%       36%       37%       38%       39%       40%       41%       42%       43%       44%       45%       46%       47%       48%       49%       50% 
#1.00      1.00      1.00      1.00      2.00      2.00      2.50      3.00      4.00      4.00      5.00      6.00      7.00      8.50     10.00     11.50     13.00 
#51%       52%       53%       54%       55%       56%       57%       58%       59%       60%       61%       62%       63%       64%       65%       66%       67% 
#15.00     17.50     20.00     23.00     26.50     30.00     34.00     38.00     42.50     48.50     54.00     60.50     67.00     74.00     82.50     91.00    101.00 
#68%       69%       70%       71%       72%       73%       74%       75%       76%       77%       78%       79%       80%       81%       82%       83%       84% 
#111.00    123.00    133.50    146.00    159.50    173.00    188.50    205.25    223.00    242.00    263.00    285.50    307.50    333.50    359.50    395.00    430.50 
#85%       86%       87%       88%       89%       90%       91%       92%       93%       94%       95%       96%       97%       98%       99%      100% 
#  469.00    509.00    558.50    604.50    664.00    730.50    810.00    899.50   1015.09   1157.22   1340.70   1619.50   2039.72   2848.50   4861.94 333303.00 

##delete the genes with high and low expression

nrow(subset)


head(design)
design$groups.vector

#fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 30)
library(MASS)
fit <- p.vector(sub.abiotic, design,family=negative.binomial(5),counts=TRUE,epsilon=0.00001)
tstep <- T.fit(fit,family= negative.binomial(5),epsilon=0.00001)

#get significant genes
sigs <- get.siggenes(tstep, vars = "all")
nrow(sigs$summary)
cluster <- see.genes(sigs$sig.genes, k = 7)
quartz.save(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/test1.pdf", bg="white",type = "pdf", device = dev.cur())
write.table(cluster$cut, file = "Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/output_cluster.txt")
write.table(sigs$sig.genes$sig.profiles, file = "Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/output_sig.profiles")
write.table(sigs$sig.genes$sig.pvalues, file = "Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/output_sig.pvalues")
sigs$summary
suma2Venn(sigs$summary[, c(2:4)])
suma2Venn(sigs$summary[, c(1:4)])

head(sigs)



STMDE66 <- data.abiotic[rownames(data.abiotic)=="Brasy1G000100.v1.1", ]
STMDE66
par(mar = c(4, 4, 4, 4))
PlotGroups (STMDE66, edesign = edesign.abiotic,xlab = "Time", ylab = "Expression value",sub="Brasy1G000100.v1.1",cex.xaxis = 1.5)
STMDE66 <- data.abiotic[rownames(data.abiotic)=="Brasy1G000200.v1.1", ]

PlotGroups (STMDE66, edesign = edesign.abiotic,xlab = "Time", ylab = "Expression value",sub="Brasy1G000200.v1.1",cex.xaxis = 1.5)

STMDE66 <- data.abiotic[rownames(data.abiotic)=="Brasy1G000200.v1.1", ]

PlotGroups (STMDE66, edesign = edesign.abiotic,xlab = "Time", ylab = "Expression value",sub="Brasy1G000200.v1.1",cex.xaxis = 1.5)



PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T, dis = design$dis, groups.vector = design$groups.vector)


names(sigs)
names(sigs$sig.genes)
names(sigs$sig.genes$ColdvsControl)
suma2Venn(sigs$summary[, c(2:4)])
suma2Venn(sigs$summary[, c(1:4)])
see.genes(sigs$sig.genes$ColdvsControl, show.fit = T, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 10)
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)
PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T, dis = design$dis, groups.vector = design$groups.vector)


################################################
edesign.abiotic <- data.frame(read.delim(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/design.abiotic",sep="\t", header=T))
rownames(edesign.abiotic) <- edesign.abiotic$X
edesign.abiotic <- edesign.abiotic[,(2:7)]
edesign.abiotic 

data <- read.delim(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/Gene_DGE_counts_time.txt", sep="\t", header=T)
data.abiotic <- data.frame(data[,c(1:61)])
rownames(data.abiotic) <- data.abiotic$GeneID
data.abiotic <- data.abiotic[,c(2:61)]
head(data.abiotic)
data.abiotic$median <- apply(data.abiotic, 1, median) 
subset <- data.abiotic[(data.abiotic$median<=1923.475 & data.abiotic$median>=10),]
head(subset)
sub.abiotic <-subset[,c(1:60)]

head(sub.abiotic)
design <- make.design.matrix(edesign.abiotic, degree = 4)

#nrow(data.abiotic)
#make the first column as row name

edesign.abiotic

quantile(data.abiotic$median,probs = seq(0,1,0.01))
#0%         1%         2%         3%         4%         5%         6%         7%         8%         9%        10%        11%        12%        13%        14% 
#0.000      0.000      0.000      0.000      0.520      1.000      1.000      1.500      2.000      2.000      2.500      3.000      3.500      4.000      4.500 
#15%        16%        17%        18%        19%        20%        21%        22%        23%        24%        25%        26%        27%        28%        29% 
#5.000      6.000      6.500      7.000      8.000      9.000     10.000     11.000     12.000     13.500     15.000     16.000     18.000     19.500     21.500 
#30%        31%        32%        33%        34%        35%        36%        37%        38%        39%        40%        41%        42%        43%        44% 
#23.500     26.000     28.000     30.500     33.000     36.000     38.500     41.500     45.000     48.500     52.000     56.000     60.500     64.500     69.500 
#45%        46%        47%        48%        49%        50%        51%        52%        53%        54%        55%        56%        57%        58%        59% 
#73.725     79.500     85.000     90.500     97.000    103.500    110.000    117.760    125.000    132.000    139.500    148.500    157.000    166.000    174.500 
#60%        61%        62%        63%        64%        65%        66%        67%        68%        69%        70%        71%        72%        73%        74% 
#185.000    195.500    206.500    218.500    230.500    243.000    256.500    270.000    285.500    299.000    314.000    332.855    349.000    369.500    392.000 
#75%        76%        77%        78%        79%        80%        81%        82%        83%        84%        85%        86%        87%        88%        89% 
#416.000    438.500    464.500    491.000    518.500    551.500    581.905    613.000    654.000    695.000    742.500    793.000    845.935    912.000    989.335 
#90%        91%        92%        93%        94%        95%        96%        97%        98%        99%       100% 
#1065.950   1170.000   1291.500   1441.000   1642.940   1923.475   2284.500   2886.745   3933.330   6772.255 333303.000 

##delete the genes with high and low expression

nrow(subset)


head(design)
design$groups.vector

#fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 30)
library(MASS)
fit <- p.vector(sub.abiotic, design,family=negative.binomial(5),counts=TRUE,epsilon=0.00001)
tstep <- T.fit(fit,family= negative.binomial(5),epsilon=0.00001)

#get significant genes
sigs <- get.siggenes(tstep, vars = "all")

cluster <- see.genes(sigs$sig.genes, k = 7)
quartz.save(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/DGE_time.pdf", bg="white",type = "pdf", device = dev.cur())
quartz.save(file="Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/DGE_treat.pdf", bg="white",type = "pdf", device = dev.cur())

write.table(cluster$cut, file = "Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/output_cluster.txt")
write.table(sigs$sig.genes$sig.profiles, file = "Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/output_sig.profiles")
write.table(sigs$sig.genes$sig.pvalues, file = "Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/output_sig.pvalues")
