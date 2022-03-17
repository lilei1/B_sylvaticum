#Written by Li Lei 2020-11-02
#This script is for calling the DGE for the cold datasets and 
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
#biocLite("httr", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
biocLite("Hmisc")
biocLite("DESeq2")

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
metaData <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/cold_tissue_meta.txt",header = T)
head(metaData)
nrow(metaData)
#9
#extract the column I needed!!!
col.num <- which( colnames(counts) %in% metaData$code)
row.num <- which( metaData$code %in% colnames(counts) )
counts <- counts[, col.num]
mata <- metaData[row.num,]
ncol(counts)
head(counts)
colnames(counts) <- metaData$exp
head(counts)
metaData$group <- paste0(metaData$tissue, ".", metaData$treatment)
#Find the best cutoff for the downstream analysis
#This is to convert a dataframe into single column.
counts_long <- gather(counts,code,counts)
head(counts_long)
#nrow(counts_long)
#this is to find the reasonable cutoff for the threshold setting
quantile(counts_long$counts,probs=seq(0,1,0.01))
#0%        1%        2%        3%        4%        5%        6%        7%        8% 
#0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00 
#9%       10%       11%       12%       13%       14%       15%       16%       17% 
#0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00 
#18%       19%       20%       21%       22%       23%       24%       25%       26% 
#0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00      0.00 
#27%       28%       29%       30%       31%       32%       33%       34%       35% 
#0.00      0.00      0.00      0.00      1.00      1.00      1.00      1.00      2.00 
#36%       37%       38%       39%       40%       41%       42%       43%       44% 
#2.00      2.00      3.00      3.00      4.00      4.00      5.00      6.00      7.00 
#45%       46%       47%       48%       49%       50%       51%       52%       53% 
#8.00      9.00     10.00     12.00     13.00     15.00     17.00     19.00     22.00 
#54%       55%       56%       57%       58%       59%       60%       61%       62% 
#24.00     27.00     30.00     34.00     37.00     41.00     46.00     50.00     55.00 
#63%       64%       65%       66%       67%       68%       69%       70%       71% 
#61.00     67.00     73.00     80.00     87.00     95.00    103.00    112.00    122.00 
#72%       73%       74%       75%       76%       77%       78%       79%       80% 
#132.00    143.00    155.00    167.00    180.00    194.00    209.00    225.00    242.00 
#81%       82%       83%       84%       85%       86%       87%       88%       89% 
#261.00    282.00    304.00    329.00    356.00    386.00    420.00    457.00    499.00 
#90%       91%       92%       93%       94%       95%       96%       97%       98% 
#549.00    606.00    675.00    760.00    867.00   1006.00   1194.00   1482.00   1976.00 
#99%      100% 
#  3203.15 459040.00 
bioreplicates.cutoff = 3# 62% of the data were above 5
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
#?rowSums
nrow(counts)
#[1] 29409
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
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 
#14%          15%          16%          17%          18%          19%          20% 
#1.013985e-01 1.051746e-01 1.138477e-01 1.203020e-01 1.226214e-01 2.027970e-01 2.183556e-01 
#21%          22%          23%          24%          25%          26%          27% 
#2.406041e-01 2.780402e-01 3.275333e-01 3.645403e-01 4.206985e-01 4.827413e-01 5.397399e-01 
#28%          29%          30%          31%          32%          33%          34% 
#6.083911e-01 6.908305e-01 7.671855e-01 8.505941e-01 9.624163e-01 1.082718e+00 1.203020e+00 
#35%          36%          37%          38%          39%          40%          41% 
#1.327539e+00 1.472445e+00 1.619220e+00 1.799133e+00 1.973663e+00 2.172336e+00 2.401911e+00 
#42%          43%          44%          45%          46%          47%          48% 
#2.618497e+00 2.849546e+00 3.127853e+00 3.397536e+00 3.696290e+00 4.019132e+00 4.355963e+00 
#49%          50%          51%          52%          53%          54%          55% 
#4.712711e+00 5.097404e+00 5.517963e+00 5.937139e+00 6.396322e+00 6.901082e+00 7.433634e+00 
#56%          57%          58%          59%          60%          61%          62% 
#7.993335e+00 8.568658e+00 9.201443e+00 9.896196e+00 1.058658e+01 1.130075e+01 1.208953e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.288043e+01 1.373102e+01 1.466293e+01 1.561242e+01 1.661241e+01 1.765748e+01 1.876108e+01 
#70%          71%          72%          73%          74%          75%          76% 
#1.991308e+01 2.114334e+01 2.245199e+01 2.380076e+01 2.525875e+01 2.674856e+01 2.839715e+01 
#77%          78%          79%          80%          81%          82%          83% 
#3.016487e+01 3.207202e+01 3.405084e+01 3.622365e+01 3.858332e+01 4.112375e+01 4.391024e+01 
#84%          85%          86%          87%          88%          89%          90% 
#4.694827e+01 5.023260e+01 5.377501e+01 5.775811e+01 6.245791e+01 6.784112e+01 7.408820e+01 
#91%          92%          93%          94%          95%          96%          97% 
#8.132062e+01 9.005353e+01 1.004783e+02 1.134649e+02 1.304588e+02 1.543333e+02 1.904272e+02 
#98%          99%         100% 
#2.513673e+02 4.070024e+02 4.827936e+04 

count.cutoff = 1 #86% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#23331
nrow(counts)
#29409
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#23331
grep("Brasy1G039900", rownames(counts))#check if this CBF got filtered
#269
#it suggested that the filtering seems good!
#21794/36927=0.6318141 only keep around 63% genes

## VST instead of voom
ncol(counts)
#creat a group 
design <- model.matrix(~0 + group,metaData)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metaData, design = design)
head(counts(dds))
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
#str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#This file is for WGCNA
write.table(x = expr.vst,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/normalized_VST_final_frezze_tissue_formal.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

###check the genes with very high variance. After checking their quantile and I found that take var = 0.3 could 
#be the best choice.
CL <- expr.vst[,1:3]
head(CL)
quantile (rowVars(CL),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0004247153 0.0010006251 0.0016209992 
#7%           8%           9%          10%          11%          12%          13% 
#0.0023603304 0.0030588345 0.0036963735 0.0044109136 0.0050769061 0.0058930136 0.0066285900 
#14%          15%          16%          17%          18%          19%          20% 
#0.0073178125 0.0081754126 0.0089297792 0.0098381672 0.0106996478 0.0114528968 0.0123529591 
#21%          22%          23%          24%          25%          26%          27% 
#0.0133313690 0.0142432883 0.0151517666 0.0160748124 0.0171051514 0.0180691253 0.0191416021 
#28%          29%          30%          31%          32%          33%          34% 
#0.0201945660 0.0213308193 0.0224202669 0.0235542299 0.0244400514 0.0255977891 0.0266268605 
#35%          36%          37%          38%          39%          40%          41% 
#0.0279022721 0.0289208469 0.0302479564 0.0316751405 0.0330143105 0.0342312941 0.0355943736 
#42%          43%          44%          45%          46%          47%          48% 
#0.0372519993 0.0387240436 0.0401798767 0.0417970852 0.0434176434 0.0452453121 0.0468906075 
#49%          50%          51%          52%          53%          54%          55% 
#0.0486675907 0.0507857462 0.0528365927 0.0549596369 0.0572831577 0.0597633268 0.0618954176 
#56%          57%          58%          59%          60%          61%          62% 
#0.0649326197 0.0675232133 0.0703133283 0.0731075013 0.0762866548 0.0792178325 0.0822701890 
#63%          64%          65%          66%          67%          68%          69% 
#0.0854903611 0.0889586039 0.0926417577 0.0969164885 0.1010756288 0.1053387559 0.1093508138 
#70%          71%          72%          73%          74%          75%          76% 
#0.1141290342 0.1193034402 0.1244242849 0.1295977909 0.1365288980 0.1443834906 0.1514561664 
#77%          78%          79%          80%          81%          82%          83% 
#0.1591495215 0.1673335190 0.1768479996 0.1859209036 0.1890079659 0.1946016181 0.2027089750 
#84%          85%          86%          87%          88%          89%          90% 
#0.2075550869 0.2208521868 0.2367049742 0.2558557544 0.2772383668 0.2951539449 0.3117907791 
#91%          92%          93%          94%          95%          96%          97% 
#0.3436759975 0.3699427100 0.3965573530 0.4382788935 0.4873387441 0.5555994273 0.6445575372 
#98%          99%         100% 
#0.7732301415 1.0813821388 8.6726716217 
CL_co <- rownames(CL[rowVars(CL)<0.3,]) #after check the quantile of the datasets and found 0.3 could be a good cutoff for the to 10% of the genes 
#with the highest varaince of the vst normalized counts.
CC <- expr.vst[,4:6]
head(CC)
quantile (rowVars(CC),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 6.100097e-05 1.972925e-04 3.160532e-04 4.595332e-04 5.964585e-04 7.366650e-04 
#7%           8%           9%          10%          11%          12%          13% 
#8.940162e-04 1.060998e-03 1.226567e-03 1.389436e-03 1.557894e-03 1.731026e-03 1.924205e-03 
#14%          15%          16%          17%          18%          19%          20% 
#2.116072e-03 2.307053e-03 2.511365e-03 2.713234e-03 2.917750e-03 3.141161e-03 3.388039e-03 
#21%          22%          23%          24%          25%          26%          27% 
#3.610840e-03 3.882336e-03 4.112417e-03 4.368142e-03 4.632337e-03 4.877680e-03 5.155826e-03 
#28%          29%          30%          31%          32%          33%          34% 
#5.439392e-03 5.709430e-03 6.057445e-03 6.379940e-03 6.672799e-03 7.046106e-03 7.389687e-03 
#35%          36%          37%          38%          39%          40%          41% 
#7.770741e-03 8.163068e-03 8.547804e-03 8.977874e-03 9.421025e-03 9.932679e-03 1.040060e-02 
#42%          43%          44%          45%          46%          47%          48% 
#1.088221e-02 1.133995e-02 1.185709e-02 1.232672e-02 1.290263e-02 1.345596e-02 1.408219e-02 
#49%          50%          51%          52%          53%          54%          55% 
#1.466161e-02 1.535472e-02 1.601558e-02 1.670631e-02 1.744833e-02 1.825977e-02 1.902795e-02 
#56%          57%          58%          59%          60%          61%          62% 
#1.988335e-02 2.077872e-02 2.161153e-02 2.262800e-02 2.353418e-02 2.451933e-02 2.561046e-02 
#63%          64%          65%          66%          67%          68%          69% 
#2.675736e-02 2.801499e-02 2.932893e-02 3.069657e-02 3.212725e-02 3.376805e-02 3.546469e-02 
#70%          71%          72%          73%          74%          75%          76% 
#3.735353e-02 3.899851e-02 4.078120e-02 4.280485e-02 4.507800e-02 4.760678e-02 5.000866e-02 
#77%          78%          79%          80%          81%          82%          83% 
#5.283474e-02 5.594384e-02 5.908404e-02 6.236528e-02 6.663853e-02 7.101572e-02 7.560817e-02 
#84%          85%          86%          87%          88%          89%          90% 
#8.034911e-02 8.602149e-02 9.166980e-02 9.914952e-02 1.064109e-01 1.155735e-01 1.242719e-01 
#91%          92%          93%          94%          95%          96%          97% 
#1.328133e-01 1.475780e-01 1.645232e-01 1.808301e-01 2.041460e-01 2.365537e-01 2.907442e-01 
#98%          99%         100% 
#3.734699e-01 5.168004e-01 4.140782e+00 
CC_co <- rownames(CC[rowVars(CC)<0.3,])

FL <- expr.vst[,7:9]
head(FL)
quantile (rowVars(FL),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 0.000000e+00 8.808496e-05 6.257241e-04 1.344575e-03 2.226587e-03 
#7%           8%           9%          10%          11%          12%          13% 
#2.970503e-03 3.765909e-03 4.606236e-03 5.469585e-03 6.341746e-03 7.162565e-03 8.080251e-03 
#14%          15%          16%          17%          18%          19%          20% 
#8.854590e-03 9.799112e-03 1.074811e-02 1.162339e-02 1.267670e-02 1.368477e-02 1.464018e-02 
#21%          22%          23%          24%          25%          26%          27% 
#1.562758e-02 1.672239e-02 1.777202e-02 1.886489e-02 1.981745e-02 2.101905e-02 2.220663e-02 
#28%          29%          30%          31%          32%          33%          34% 
#2.333559e-02 2.461232e-02 2.577394e-02 2.687463e-02 2.818674e-02 2.939905e-02 3.048845e-02 
#35%          36%          37%          38%          39%          40%          41% 
#3.176826e-02 3.285311e-02 3.428487e-02 3.577342e-02 3.722742e-02 3.840201e-02 3.994541e-02 
#42%          43%          44%          45%          46%          47%          48% 
#4.177812e-02 4.346158e-02 4.515026e-02 4.689750e-02 4.833854e-02 5.002133e-02 5.201557e-02 
#49%          50%          51%          52%          53%          54%          55% 
#5.383519e-02 5.621673e-02 5.817343e-02 6.032785e-02 6.249694e-02 6.462480e-02 6.684577e-02 
#56%          57%          58%          59%          60%          61%          62% 
#6.939451e-02 7.182936e-02 7.442177e-02 7.711021e-02 7.976921e-02 8.260601e-02 8.580661e-02 
#63%          64%          65%          66%          67%          68%          69% 
#8.898553e-02 9.219802e-02 9.559415e-02 9.925845e-02 1.033837e-01 1.070927e-01 1.118110e-01 
#70%          71%          72%          73%          74%          75%          76% 
#1.164000e-01 1.205448e-01 1.256396e-01 1.309222e-01 1.360191e-01 1.429100e-01 1.492192e-01 
#77%          78%          79%          80%          81%          82%          83% 
#1.569756e-01 1.639774e-01 1.719936e-01 1.783507e-01 1.815304e-01 1.836998e-01 1.876154e-01 
#84%          85%          86%          87%          88%          89%          90% 
#1.960937e-01 2.094803e-01 2.223304e-01 2.368126e-01 2.521738e-01 2.738205e-01 2.870185e-01 
#91%          92%          93%          94%          95%          96%          97% 
#3.051547e-01 3.327491e-01 3.584912e-01 3.856371e-01 4.239594e-01 4.829746e-01 5.429037e-01 
#98%          99%         100% 
#6.457258e-01 8.801205e-01 1.424281e+01 
FL_co <- rownames(FL[rowVars(FL)<0.3,])
#at least 1 treatment has variance <0.01
FC <- expr.vst[,10:12]
head(FC)
quantile (rowVars(FC),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.0000000000 0.0001238319 0.0003219687 0.0005277496 0.0007336517 0.0009322888 0.0011658269 
#7%           8%           9%          10%          11%          12%          13% 
#0.0013860280 0.0016271666 0.0018374792 0.0020519157 0.0022995876 0.0025391238 0.0027896254 
#14%          15%          16%          17%          18%          19%          20% 
#0.0030877888 0.0033446606 0.0035958746 0.0038810054 0.0041508560 0.0044985883 0.0047870486 
#21%          22%          23%          24%          25%          26%          27% 
#0.0051079765 0.0054141888 0.0057654796 0.0061076470 0.0064868108 0.0068434398 0.0072220534 
#28%          29%          30%          31%          32%          33%          34% 
#0.0076016848 0.0079903342 0.0084009331 0.0088157715 0.0092405988 0.0097338389 0.0102009588 
#35%          36%          37%          38%          39%          40%          41% 
#0.0106095282 0.0110438792 0.0116442865 0.0121334778 0.0127161569 0.0133324214 0.0138982993 
#42%          43%          44%          45%          46%          47%          48% 
#0.0145784070 0.0152223254 0.0157815616 0.0164578982 0.0172124938 0.0178366799 0.0184922670 
#49%          50%          51%          52%          53%          54%          55% 
#0.0193767125 0.0200928150 0.0210992884 0.0219222399 0.0227794412 0.0236582574 0.0246533581 
#56%          57%          58%          59%          60%          61%          62% 
#0.0256635804 0.0266147862 0.0276524483 0.0288353734 0.0300297090 0.0314984706 0.0328594639 
#63%          64%          65%          66%          67%          68%          69% 
#0.0343644068 0.0358246255 0.0376053713 0.0393353752 0.0410530759 0.0430644623 0.0451019592 
#70%          71%          72%          73%          74%          75%          76% 
#0.0474117754 0.0495360164 0.0517591288 0.0545953107 0.0576022862 0.0604641133 0.0635217428 
#77%          78%          79%          80%          81%          82%          83% 
#0.0669904254 0.0707934551 0.0744625746 0.0783640098 0.0825278824 0.0874648521 0.0933191877 
#84%          85%          86%          87%          88%          89%          90% 
#0.0990186484 0.1045737455 0.1107510967 0.1175718092 0.1268421539 0.1361001532 0.1459002988 
#91%          92%          93%          94%          95%          96%          97% 
#0.1575949240 0.1733207833 0.1910684669 0.2134826219 0.2394256552 0.2775452400 0.3287367004 
#98%          99%         100% 
#0.4012188355 0.5628351583 8.2236554296 
FC_co <- rownames(FC[rowVars(FC)<0.3,])
RL <- expr.vst[,13:15]
head(RL)
quantile (rowVars(RL),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0005496672 0.0009676443 0.0015549859 
#7%           8%           9%          10%          11%          12%          13% 
#0.0021669847 0.0025164923 0.0030745909 0.0035652743 0.0041038671 0.0046699793 0.0052936035 
#14%          15%          16%          17%          18%          19%          20% 
#0.0058900641 0.0064929566 0.0070421093 0.0076815003 0.0082631570 0.0088389309 0.0095792066 
#21%          22%          23%          24%          25%          26%          27% 
#0.0102398860 0.0109830070 0.0116826781 0.0123928887 0.0130349923 0.0138610919 0.0146694185 
#28%          29%          30%          31%          32%          33%          34% 
#0.0155094812 0.0162886968 0.0172037297 0.0181385549 0.0190358282 0.0200116808 0.0210701264 
#35%          36%          37%          38%          39%          40%          41% 
#0.0220911849 0.0229161198 0.0240051743 0.0251263294 0.0262001470 0.0274137479 0.0284538038 
#42%          43%          44%          45%          46%          47%          48% 
#0.0296623268 0.0307809354 0.0322076979 0.0336592125 0.0350177863 0.0365732707 0.0381321980 
#49%          50%          51%          52%          53%          54%          55% 
#0.0396623062 0.0410912570 0.0427906339 0.0443756447 0.0461131817 0.0480136092 0.0501206773 
#56%          57%          58%          59%          60%          61%          62% 
#0.0521334225 0.0543796508 0.0564541001 0.0589301729 0.0614018138 0.0643112203 0.0675667725 
#63%          64%          65%          66%          67%          68%          69% 
#0.0704747591 0.0734602652 0.0767575426 0.0798650473 0.0830902351 0.0865433573 0.0906852028 
#70%          71%          72%          73%          74%          75%          76% 
#0.0952365246 0.0990287513 0.1039473081 0.1089396729 0.1139866505 0.1202133798 0.1265239903 
#77%          78%          79%          80%          81%          82%          83% 
#0.1335393142 0.1411120049 0.1488604129 0.1577830802 0.1647685276 0.1715382983 0.1817151340 
#84%          85%          86%          87%          88%          89%          90% 
#0.1893114545 0.1984582516 0.2054615980 0.2109996119 0.2244830559 0.2446610798 0.2670467810 
#91%          92%          93%          94%          95%          96%          97% 
#0.2875226984 0.3162260283 0.3367688608 0.3703623087 0.4120062849 0.4710726096 0.5439769622 
#98%          99%         100% 
#0.6551190489 0.9299135923 9.2313054425
RL_co <- rownames(RL[rowVars(RL)<0.3,])

RC <- expr.vst[,16:18]
head(RC)
quantile (rowVars(RC),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.0000000000 0.0004820726 0.0010190321 0.0016039526 0.0021642203 0.0027592490 0.0033992549 
#7%           8%           9%          10%          11%          12%          13% 
#0.0040274466 0.0046917746 0.0053793320 0.0060683647 0.0067815941 0.0075461301 0.0082575924 
#14%          15%          16%          17%          18%          19%          20% 
#0.0089590867 0.0096780027 0.0104058631 0.0111942042 0.0120531677 0.0128435706 0.0137090282 
#21%          22%          23%          24%          25%          26%          27% 
#0.0145715043 0.0154898061 0.0164390452 0.0173981726 0.0183684444 0.0194485219 0.0204075592 
#28%          29%          30%          31%          32%          33%          34% 
#0.0214830906 0.0226552826 0.0237925868 0.0249778158 0.0259565497 0.0270790480 0.0284168436 
#35%          36%          37%          38%          39%          40%          41% 
#0.0296887520 0.0308142380 0.0320964039 0.0335567632 0.0349261051 0.0361229506 0.0375107186 
#42%          43%          44%          45%          46%          47%          48% 
#0.0388034981 0.0404434700 0.0420174501 0.0435149925 0.0451461708 0.0469117110 0.0486157218 
#49%          50%          51%          52%          53%          54%          55% 
#0.0502943686 0.0519819452 0.0537685070 0.0558658988 0.0578407155 0.0597614907 0.0618141147 
#56%          57%          58%          59%          60%          61%          62% 
#0.0639533906 0.0660335384 0.0682516360 0.0706502930 0.0727212010 0.0753197305 0.0779305216 
#63%          64%          65%          66%          67%          68%          69% 
#0.0804674696 0.0828757290 0.0856669669 0.0884012194 0.0914164497 0.0948316053 0.0991079470 
#70%          71%          72%          73%          74%          75%          76% 
#0.1021786039 0.1058821655 0.1098893891 0.1143777949 0.1186654895 0.1228139390 0.1268616430 
#77%          78%          79%          80%          81%          82%          83% 
#0.1313619434 0.1362364906 0.1420128011 0.1464782486 0.1527172384 0.1594245104 0.1652044290 
#84%          85%          86%          87%          88%          89%          90% 
#0.1729314045 0.1811468614 0.1905991063 0.2001982402 0.2105711439 0.2231691900 0.2365844026 
#91%          92%          93%          94%          95%          96%          97% 
#0.2508386440 0.2714798305 0.2929986827 0.3223181322 0.3561548158 0.3969075575 0.4596484017 
#98%          99%         100% 
#0.5565251815 0.7593637280 5.1174295192 
RC_co <- rownames(RC[rowVars(RC)<0.3,])

low_var_gene <- union(union(union(union(union(CC_co,CL_co),FC_co),FL_co),RC_co),RL_co)# here are the genes with loc varaince for each treatment 
length(low_var_gene)
#23310
counts <- counts[rownames(counts) %in% low_var_gene,]
nrow(counts)
#23310
grep("Brasy1G039900", rownames(counts))#check if this CBF got filtered
#269
#it suggested that the filtering seems good!
#23310/36927=0.6312454 only keep around 60% genes

## VST instead of voom
#creat a group 
design <- model.matrix(~0 + group,metaData)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = metaData, design = design)
head(counts(dds))
#colSums(counts(dds)) %>% barplot
expr.vst <- assay(DESeq2::vst(dds))
head(expr.vst)
#str(expr.vst)
expr.vst <- as.data.frame(expr.vst)
#This file is for WGCNA
write.table(x = expr.vst,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/normalized_VST_final_freeze_tissue_novargen.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

####
ddsTC <- DESeq(dds)
#heatmap:
library("pheatmap")
rld <- vst( ddsTC )
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/PCA_freeze_tissue_syl_novargen.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
plotPCA(rld, intgroup="group")
dev.off()
colData(dds)
head(rld)
library( "genefilter" )
#install.packages("gplots")
library(gplots)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 1000 )
head(topVarGenes)
#dev.off()
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_1000_DGE_freeze_syl_novargen.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( cold="#11C5F6", freeze="#1149F6", recovery="gray")[
             colData(rld)$treatment ] )
dev.off()

##cold versus freeze leaf
colData(ddsTC)
res_freeze_leaf <- results(ddsTC, contrast=list("groupleaf.freeze", "groupleaf.cold"))
EnhancedVolcano(res_freeze_leaf,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'freeze_leaf',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/freeze_leaf_syl.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#recovery vs freeze leaf
res_recovery_leaf <- results(ddsTC, contrast=list("groupleaf.recovery", "groupleaf.freeze"))
EnhancedVolcano(res_recovery_leaf,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'recovery_leaf',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/recovery_leaf_syl.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#recovery versus cold leaf
res_recoverycold_leaf <- results(ddsTC, contrast=list("groupleaf.recovery", "groupleaf.cold"))
EnhancedVolcano(res_recoverycold_leaf,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'recovery vs cold leaf',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/recovery_vs_cold_leaf_syl.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
##crown
##cold versus freeze crown
colData(ddsTC)
res_freeze_crown <- results(ddsTC, contrast=list("groupcrown.freeze", "groupcrown.cold"))
EnhancedVolcano(res_freeze_crown,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'freeze_crown',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/freeze_crown_syl.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#recovery vs freeze crown
res_recovery_crown <- results(ddsTC, contrast=list("groupcrown.recovery", "groupcrown.freeze"))
EnhancedVolcano(res_recovery_crown,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'recovery_crown',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/recovery_crown_syl.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#recovery versus cold crown
res_recoverycold_crown <- results(ddsTC, contrast=list("groupcrown.recovery", "groupcrown.cold"))
EnhancedVolcano(res_recoverycold_crown,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'recovery vs cold crown',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=1,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/recovery_cold_crown_syl.pdf",
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
#cold
res_freeze_leaf <- res_freeze_leaf[!is.na(res_freeze_leaf$padj),]
sigres_res_freeze_leaf_up <- res_freeze_leaf[(res_freeze_leaf$padj<0.05 & res_freeze_leaf$log2FoldChange>= 1),]
nb_res_freeze_leaf_up <- nrow(sigres_res_freeze_leaf_up)
#693
write.table(x = sigres_res_freeze_leaf_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_freeze_leaf_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_res_freeze_leaf_down <- res_freeze_leaf[(res_freeze_leaf$padj<0.05 & res_freeze_leaf$log2FoldChange<= -1),]
nb_res_freeze_leaf_down <- nrow(sigres_res_freeze_leaf_down)
#119
write.table(x = sigres_res_freeze_leaf_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_freeze_leaf_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_res_freeze_leaf_DEG <- res_freeze_leaf[(res_freeze_leaf$padj<0.05 & abs(res_freeze_leaf$log2FoldChange) >= 1 ),]
nb_sigres_res_freeze_leaf_DEG <- nrow(sigres_res_freeze_leaf_DEG)
#812
#set the logfoldchanges as 5 and extract those genes
sigres_res_freeze_leaf_DEG_strict <- res_freeze_leaf[(res_freeze_leaf$padj<0.05 & abs(res_freeze_leaf$log2FoldChange) >= 2),]
head(sigres_res_freeze_leaf_DEG_strict)
nrow(sigres_res_freeze_leaf_DEG_strict)
#475
write.table(x = sigres_res_freeze_leaf_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_freeze_leaf_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

##recovery from freeze
res_recovery_leaf <- res_recovery_leaf[!is.na(res_recovery_leaf$padj),]
sigres_res_recovery_leaf_up <- res_recovery_leaf[(res_recovery_leaf$padj<0.05 & res_recovery_leaf$log2FoldChange>= 1),]
nb_res_recovery_leaf_up <- nrow(sigres_res_recovery_leaf_up)
#2123
write.table(x = sigres_res_recovery_leaf_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recovery_leaf_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_res_recovery_leaf_down <- res_recovery_leaf[(res_recovery_leaf$padj<0.05 & res_recovery_leaf$log2FoldChange <= -1),]
nb_res_recovery_leaf_down <- nrow(sigres_res_recovery_leaf_down)
#2345
write.table(x = sigres_res_recovery_leaf_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recovery_leaf_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_res_recovery_leaf_DEG <- res_recovery_leaf[(res_recovery_leaf$padj<0.05 & abs(res_recovery_leaf$log2FoldChange) >= 1 ),]
nb_sigres_res_recovery_leaf_DEG <- nrow(sigres_res_recovery_leaf_DEG)
#4468

#set the logfoldchanges as 5 and extract those genes
sigres_res_recovery_leaf_DEG_strict <- res_recovery_leaf[(res_recovery_leaf$padj<0.05 & abs(res_recovery_leaf$log2FoldChange) >= 2),]
head(sigres_res_recovery_leaf_DEG_strict)
nrow(sigres_res_recovery_leaf_DEG_strict)
#1775
write.table(x = sigres_res_recovery_leaf_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recovery_leaf_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#recovery versus the cold
res_recoverycold_leaf <- res_recoverycold_leaf[!is.na(res_recoverycold_leaf$padj),]
sigres_res_recoverycold_leaf_up <- res_recoverycold_leaf[(res_recoverycold_leaf$padj<0.05 & res_recoverycold_leaf$log2FoldChange>= 1),]
nb_res_recoverycold_leaf_up <- nrow(sigres_res_recoverycold_leaf_up)
#2437
write.table(x = sigres_res_recoverycold_leaf_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recoverycold_leaf_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_res_recoverycold_leaf_down <- res_recoverycold_leaf[(res_recoverycold_leaf$padj<0.05 & res_recoverycold_leaf$log2FoldChange<= -1),]
nb_res_recoverycold_leaf_down <- nrow(sigres_res_recoverycold_leaf_down)
#2346
write.table(x = sigres_res_recoverycold_leaf_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recoverycold_leaf_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_res_recoverycold_leaf_DEG <- res_recoverycold_leaf[(res_recoverycold_leaf$padj<0.05 & abs(res_recoverycold_leaf$log2FoldChange) >= 1 ),]
nb_sigres_res_recoverycold_leaf_DEG <- nrow(sigres_res_recoverycold_leaf_DEG)
#4783
#set the logfoldchanges as 5 and extract those genes
sigres_res_recoverycold_leaf_DEG_strict <- res_recoverycold_leaf[(res_recoverycold_leaf$padj<0.05 & abs(res_recoverycold_leaf$log2FoldChange) >= 2),]
head(sigres_res_recoverycold_leaf_DEG_strict)
nrow(sigres_res_recoverycold_leaf_DEG_strict)
#2068
write.table(x = sigres_res_recoverycold_leaf_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recoverycold_leaf_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
###Crown
res_freeze_crown <- res_freeze_crown[!is.na(res_freeze_crown$padj),]
sigres_res_freeze_crown_up <- res_freeze_crown[(res_freeze_crown$padj<0.05 & res_freeze_crown$log2FoldChange>= 1),]
nb_res_freeze_crown_up <- nrow(sigres_res_freeze_crown_up)
#1246
write.table(x = sigres_res_freeze_crown_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_freeze_crown_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_res_freeze_crown_down <- res_freeze_crown[(res_freeze_crown$padj<0.05 & res_freeze_crown$log2FoldChange<= -1),]
nb_res_freeze_crown_down <- nrow(sigres_res_freeze_crown_down)
#586
write.table(x = sigres_res_freeze_crown_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_freeze_crown_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_res_freeze_crown_DEG <- res_freeze_crown[(res_freeze_crown$padj<0.05 & abs(res_freeze_crown$log2FoldChange) >= 1 ),]
nb_sigres_res_freeze_crown_DEG <- nrow(sigres_res_freeze_crown_DEG)
#812
#set the logfoldchanges as 5 and extract those genes
sigres_res_freeze_crown_DEG_strict <- res_freeze_crown[(res_freeze_crown$padj<0.05 & abs(res_freeze_crown$log2FoldChange) >= 2),]
head(sigres_res_freeze_crown_DEG_strict)
nrow(sigres_res_freeze_crown_DEG_strict)
#838
write.table(x = sigres_res_freeze_crown_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_freeze_crown_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

##recovery from freeze
res_recovery_crown <- res_recovery_crown[!is.na(res_recovery_crown$padj),]
sigres_res_recovery_crown_up <- res_recovery_crown[(res_recovery_crown$padj<0.05 & res_recovery_crown$log2FoldChange>= 1),]
nb_res_recovery_crown_up <- nrow(sigres_res_recovery_crown_up)
#2542
write.table(x = sigres_res_recovery_crown_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recovery_crown_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_res_recovery_crown_down <- res_recovery_crown[(res_recovery_crown$padj<0.05 & res_recovery_crown$log2FoldChange <= -1),]
nb_res_recovery_crown_down <- nrow(sigres_res_recovery_crown_down)
#1355
write.table(x = sigres_res_recovery_crown_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recovery_crown_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_res_recovery_crown_DEG <- res_recovery_crown[(res_recovery_crown$padj<0.05 & abs(res_recovery_crown$log2FoldChange) >= 1 ),]
nb_sigres_res_recovery_crown_DEG <- nrow(sigres_res_recovery_crown_DEG)
#3897

#set the logfoldchanges as 5 and extract those genes
sigres_res_recovery_crown_DEG_strict <- res_recovery_crown[(res_recovery_crown$padj<0.05 & abs(res_recovery_crown$log2FoldChange) >= 2),]
head(sigres_res_recovery_crown_DEG_strict)
nrow(sigres_res_recovery_crown_DEG_strict)
#1379
write.table(x = sigres_res_recovery_crown_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recovery_crown_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#recovery versus the cold
res_recoverycold_crown <- res_recoverycold_crown[!is.na(res_recoverycold_crown$padj),]
sigres_res_recoverycold_crown_up <- res_recoverycold_crown[(res_recoverycold_crown$padj<0.05 & res_recoverycold_crown$log2FoldChange>= 1),]
nb_res_recoverycold_crown_up <- nrow(sigres_res_recoverycold_crown_up)
#3502
write.table(x = sigres_res_recoverycold_crown_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recoverycold_crown_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
sigres_res_recoverycold_crown_down <- res_recoverycold_crown[(res_recoverycold_crown$padj<0.05 & res_recoverycold_crown$log2FoldChange<= -1),]
nb_res_recoverycold_crown_down <- nrow(sigres_res_recoverycold_crown_down)
#2907
write.table(x = sigres_res_recoverycold_crown_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recoverycold_crown_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


sigres_res_recoverycold_crown_DEG <- res_recoverycold_crown[(res_recoverycold_crown$padj<0.05 & abs(res_recoverycold_crown$log2FoldChange) >= 1 ),]
nb_sigres_res_recoverycold_crown_DEG <- nrow(sigres_res_recoverycold_crown_DEG)
#6409
#set the logfoldchanges as 5 and extract those genes
sigres_res_recoverycold_crown_DEG_strict <- res_recoverycold_crown[(res_recoverycold_crown$padj<0.05 & abs(res_recoverycold_crown$log2FoldChange) >= 2),]
head(sigres_res_recoverycold_crown_DEG_strict)
nrow(sigres_res_recoverycold_crown_DEG_strict)
#2763
write.table(x = sigres_res_recoverycold_crown_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_recoverycold_crown_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#top gene list:
top_gene_list_str <- union(union(union(union(union(rownames(sigres_res_freeze_leaf_DEG_strict),rownames(sigres_res_recovery_leaf_DEG_strict)),
                                             rownames(sigres_res_recoverycold_leaf_DEG_strict)),
                                       rownames(sigres_res_freeze_crown_DEG_strict)),
                                 rownames(sigres_res_recovery_crown_DEG_strict)),
                           rownames(sigres_res_recoverycold_crown_DEG_strict))

length(top_gene_list_str)
#4297
top_gene_list <- union(union(union(union(union(rownames(sigres_res_freeze_leaf_DEG),rownames(sigres_res_recovery_leaf_DEG)),
                                             rownames(sigres_res_recoverycold_leaf_DEG)),
                                       rownames(sigres_res_freeze_crown_DEG)),
                                 rownames(sigres_res_recovery_crown_DEG)),
                           rownames(sigres_res_recoverycold_crown_DEG))

length(top_gene_list)
#9475

#plot the genes with logfolderchanges >1
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_4297_DGE_freeze_syl_novargen#.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)
head(assay(rld))
reorder_rld <- assay(rld)[,c(1,2,3,7,8,9,13,14,15,4,5,6,10,11,12,16,17,18)]
head(reorder_rld)
heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( leaf ="green", crown="gold")[
             colData(rld)$tissue ] )
heatmap.2( assay(rld)[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( leaf ="green", crown="gold")[
             colData(rld)$tissue ] )
dev.off()
#plot the genes with logfolderchanges >2
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_9475_DGE_freeze_syl_novargen.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( reorder_rld[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( cold="#1149F6", freeze="#11C5F6", recovery="gray")[
             colData(rld)$treatment ] )
dev.off()

###rlog normalization
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_110_strictDGE_cold_syl_rlog.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rlog_ds)[row.names(assay(rlog_ds)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", cold="#11C5F6", freeze="#1149F6")[
             colData(rld)$treatment ] )
dev.off()
#rlog_ds from the plot, rlog is not great for normalization!!!!!!


pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_110_strictDGE_cold_syl_logcount.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
##Test another way
#selected <- rownames(sig)
output <- counts(dds,normalized=TRUE)[rownames(dds) %in% top_gene_list,]
output <- log2(output)
head(output)
#output[is.na(output)] = 0
#output_no_NA <- output[!is.na(resOrdered$padj),]
#na.dist <- function(x) {
#  t.dist <- dist(x)
#  t.dist <- as.matrix(t.dist)
#  t.limit <- 1.1*max(t.dist,na.rm=T)
#  t.dist[is.na(t.dist)] <- t.limit
##  t.dist <- as.dist(t.dist)
# return(t.dist)
#}
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
output[!is.finite(output)] <- 0

palette <- colorRampPalette(c("red","white","blue"))(256)
#docs <- dist( as.matrix(output), method = "euclidean")
#hclust_dist<- as.dist(docs)
#hclust_dist[is.na(hclust_dist)] <- 0
#hclust_dist[is.nan(hclust_dist)] <- 0
#sum(is.infinite(hclust_dist))  # THIS SHOULD BE 0
#h <- hclust(hclust_dist, "ward.D2")

heatmap.2(output,Colv=FALSE,
          distfun=dist_no_na,
          col = palette, scale="row", 
          dendrogram="row",
          trace="none",margin=c(4,6), 
          cexRow=0.5, cexCol=1, keysize=1 )
dev.off()
#check
Brasy8G092900.v1.1
counts[(row.names(counts) == "Brasy8G092900.v1.1"),]
#head(count)
assay(rlog_ds)[row.names(assay(rlog_ds)) == "Brasy8G092900.v1.1",]
