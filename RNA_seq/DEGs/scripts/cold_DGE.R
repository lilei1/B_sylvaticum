Written by Li Lei 2020-11-02
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
metaData <- read.delim("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/cold_meta.txt",header = T)
head(metaData)
nrow(metaData)
#18
#extract the column I needed!!!
col.num <- which( colnames(counts) %in% metaData$code)
row.num <- which( metaData$code %in% colnames(counts) )
counts <- counts[, col.num]
mata <- metaData[row.num,]
ncol(counts)
head(counts)
colnames(counts) <- metaData$exp
head(counts)
metaData$group <- paste0(metaData$acc, ".", metaData$treatment)

#Find the best cutoff for the downstream analysis
#This is to convert a dataframe into single column.
counts_long <- gather(counts,code,counts)
head(counts_long)
#nrow(counts_long)
#this is to find the reasonable cutoff for the threshold setting
quantile(counts_long$counts,probs=seq(0,1,0.01))
#0%     1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11%    12% 
#0      0      0      0      0      0      0      0      0      0      0      0      0 
#13%    14%    15%    16%    17%    18%    19%    20%    21%    22%    23%    24%    25% 
#0      0      0      0      0      0      0      0      0      0      0      0      0 
#26%    27%    28%    29%    30%    31%    32%    33%    34%    35%    36%    37%    38% 
#0      0      0      0      1      1      1      1      2      2      2      3      4 
#39%    40%    41%    42%    43%    44%    45%    46%    47%    48%    49%    50%    51% 
#5      5      7      8      9     11     13     15     17     20     22     26     29 
#52%    53%    54%    55%    56%    57%    58%    59%    60%    61%    62%    63%    64% 
#33     38     43     48     54     61     68     75     84     93    103    114    126 
#65%    66%    67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
#139    153    168    184    201    220    240    261    284    308    334    361    391 
#78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88%    89%    90% 
#423    457    494    534    578    625    677    733    797    867    944   1034   1134 
#91%    92%    93%    94%    95%    96%    97%    98%    99%   100% 
#  1256   1399   1574   1795   2074   2458   3028   4022   6314 470977 

count.cutoff = 5 # 69% of the data were above 5
bioreplicates.cutoff = 3
## RAW COUNTS ##
keep <- rowSums(counts >= count.cutoff) >= bioreplicates.cutoff
counts <- counts[keep, ]
colnames(counts)
#?rowSums
nrow(counts)
#[1] 25898
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
#5.281276e-02 6.750150e-02 9.936210e-02 1.266693e-01 1.651602e-01 2.025045e-01 2.436127e-01 
#14%          15%          16%          17%          18%          19%          20% 
#2.820123e-01 3.269296e-01 3.762148e-01 4.258345e-01 4.842851e-01 5.464916e-01 6.083349e-01 
#21%          22%          23%          24%          25%          26%          27% 
#6.777482e-01 7.600158e-01 8.450042e-01 9.450210e-01 1.035354e+00 1.144254e+00 1.253834e+00 
#28%          29%          30%          31%          32%          33%          34% 
#1.393362e+00 1.524934e+00 1.665399e+00 1.822568e+00 2.006885e+00 2.185966e+00 2.372506e+00 
#35%          36%          37%          38%          39%          40%          41% 
#2.579758e+00 2.798341e+00 3.045159e+00 3.303203e+00 3.577579e+00 3.863414e+00 4.172208e+00 
#42%          43%          44%          45%          46%          47%          48% 
#4.501478e+00 4.858774e+00 5.225719e+00 5.615446e+00 6.042283e+00 6.484104e+00 6.935018e+00 
#49%          50%          51%          52%          53%          54%          55% 
#7.430188e+00 7.960468e+00 8.509581e+00 9.091632e+00 9.690202e+00 1.035130e+01 1.102023e+01 
#56%          57%          58%          59%          60%          61%          62% 
#1.176947e+01 1.253170e+01 1.332253e+01 1.411337e+01 1.498533e+01 1.591057e+01 1.684027e+01 
#63%          64%          65%          66%          67%          68%          69% 
#1.784333e+01 1.888712e+01 1.998438e+01 2.112177e+01 2.229380e+01 2.353893e+01 2.484050e+01 
#70%          71%          72%          73%          74%          75%          76% 
#2.622058e+01 2.765004e+01 2.909785e+01 3.067690e+01 3.237422e+01 3.415840e+01 3.606188e+01 
#77%          78%          79%          80%          81%          82%          83% 
#3.807015e+01 4.021715e+01 4.247788e+01 4.500121e+01 4.768992e+01 5.053901e+01 5.359431e+01 
#84%          85%          86%          87%          88%          89%          90% 
#5.698585e+01 6.060058e+01 6.469563e+01 6.912103e+01 7.415593e+01 8.007324e+01 8.687425e+01 
#91%          92%          93%          94%          95%          96%          97% 
#9.471995e+01 1.041008e+02 1.152288e+02 1.292977e+02 1.475195e+02 1.719371e+02 2.089136e+02 
#98%          99%         100% 
#2.730608e+02 4.300059e+02 1.772088e+04
count.cutoff = 1 #90% of the data is above 1
bioreplicates.cutoff = 3
keep <- rowSums(normalized.counts >= count.cutoff) >= bioreplicates.cutoff
#abiotic_shoot_count(keep)
normalized.counts <- normalized.counts[keep, ]
nrow(normalized.counts)
#21293
nrow(counts)
#25898
#keep the datapoint pass the threshold
counts <- counts[row.names(normalized.counts),]
nrow(counts)
#22622
grep("Brasy1G039900", rownames(counts))#check if this CBF got filtered
#263
#it suggested that the filtering seems good!
#21794/36927=0.6126141 only keep around 60% genes

## VST instead of voom
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
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/normalized_VST_final_cold_formal#.txt",
            sep = "\t",
            eol = "\n",
            col.names = TRUE,
            row.names = TRUE
)

###check the genes with very high variance
Osl1_CK <- expr.vst[,1:3]
#head(CK)
quantile (rowVars(Osl1_CK),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000979966 
#7%           8%           9%          10%          11%          12%          13% 
#0.0002195554 0.0003265946 0.0004330755 0.0005679642 0.0006866593 0.0008128689 0.0009179842 
#14%          15%          16%          17%          18%          19%          20% 
#0.0010477165 0.0011598948 0.0012908153 0.0014161438 0.0015426527 0.0016705140 0.0018027352 
#21%          22%          23%          24%          25%          26%          27% 
#0.0019454866 0.0020963323 0.0022333225 0.0023999211 0.0025469281 0.0027040362 0.0028482551 
#28%          29%          30%          31%          32%          33%          34% 
#0.0029940724 0.0031536229 0.0033241660 0.0034726597 0.0036351739 0.0038149182 0.0040045260 
#35%          36%          37%          38%          39%          40%          41% 
#0.0041644272 0.0043757728 0.0045458412 0.0047313015 0.0049143768 0.0051019547 0.0053230594 
#42%          43%          44%          45%          46%          47%          48% 
#0.0055388636 0.0057335971 0.0059425394 0.0061642776 0.0063700153 0.0065939665 0.0068311460 
#49%          50%          51%          52%          53%          54%          55% 
#0.0070819304 0.0073232519 0.0074934264 0.0077590616 0.0079016589 0.0081276676 0.0083400526 
#56%          57%          58%          59%          60%          61%          62% 
#0.0086176712 0.0089075952 0.0092110375 0.0095102668 0.0097992708 0.0101123571 0.0104500783 
#63%          64%          65%          66%          67%          68%          69% 
#0.0107968519 0.0111181543 0.0115264533 0.0119054186 0.0122792775 0.0126211139 0.0130353078 
#70%          71%          72%          73%          74%          75%          76% 
#0.0134769041 0.0140495563 0.0144735704 0.0149508412 0.0154807274 0.0159847056 0.0165529156 
#77%          78%          79%          80%          81%          82%          83% 
#0.0172026159 0.0179068799 0.0187040890 0.0194984625 0.0202973786 0.0212187949 0.0221681857 
#84%          85%          86%          87%          88%          89%          90% 
#0.0231510812 0.0243422215 0.0256265055 0.0270942928 0.0285979350 0.0301422683 0.0320761385 
#91%          92%          93%          94%          95%          96%          97% 
#0.0342201633 0.0368598420 0.0401500492 0.0436055295 0.0485990659 0.0555171518 0.0654817164 
#98%          99%         100% 
#0.0812623882 0.1279247151 2.0889021260 
Osl1_CK_co <- rownames(Osl1_CK[rowVars(Osl1_CK)<0.07,]) #after check the quantile of the datasets and found 0.07 could be a good cutoff for the to 10% of the genes 
#with the highest varaince of the vst normalized counts.
Osl1_cold <- expr.vst[,4:6]
head(Osl1_cold)
quantile (rowVars(Osl1_cold),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 8.901635e-05 1.993588e-04 
#7%           8%           9%          10%          11%          12%          13% 
#3.202461e-04 4.268503e-04 5.560369e-04 6.737172e-04 7.877099e-04 9.040949e-04 1.032194e-03 
#14%          15%          16%          17%          18%          19%          20% 
#1.160678e-03 1.294187e-03 1.436995e-03 1.586973e-03 1.722519e-03 1.866355e-03 2.011875e-03 
#21%          22%          23%          24%          25%          26%          27% 
#2.151428e-03 2.305186e-03 2.432225e-03 2.615007e-03 2.767806e-03 2.945927e-03 3.122580e-03 
#28%          29%          30%          31%          32%          33%          34% 
#3.295327e-03 3.459521e-03 3.640940e-03 3.818712e-03 3.984616e-03 4.192530e-03 4.296646e-03 
#35%          36%          37%          38%          39%          40%          41% 
#4.490363e-03 4.690132e-03 4.912943e-03 5.110941e-03 5.349731e-03 5.550123e-03 5.719480e-03 
#42%          43%          44%          45%          46%          47%          48% 
#5.954026e-03 6.201905e-03 6.419940e-03 6.618882e-03 6.852852e-03 7.076980e-03 7.313459e-03 
#49%          50%          51%          52%          53%          54%          55% 
#7.565712e-03 7.694530e-03 7.915264e-03 8.183693e-03 8.457595e-03 8.765821e-03 9.086778e-03 
#56%          57%          58%          59%          60%          61%          62% 
#9.431157e-03 9.781656e-03 1.012988e-02 1.048869e-02 1.085186e-02 1.127176e-02 1.164631e-02 
#63%          64%          65%          66%          67%          68%          69% 
#1.203409e-02 1.240990e-02 1.279999e-02 1.324337e-02 1.372236e-02 1.422851e-02 1.482411e-02 
#70%          71%          72%          73%          74%          75%          76% 
#1.537487e-02 1.595753e-02 1.652669e-02 1.716921e-02 1.787517e-02 1.858338e-02 1.932719e-02 
#77%          78%          79%          80%          81%          82%          83% 
#2.014252e-02 2.104165e-02 2.206725e-02 2.304059e-02 2.401510e-02 2.515305e-02 2.645107e-02 
#84%          85%          86%          87%          88%          89%          90% 
#2.790045e-02 2.932316e-02 3.086341e-02 3.278715e-02 3.506382e-02 3.727745e-02 4.013417e-02 
#91%          92%          93%          94%          95%          96%          97% 
#4.322843e-02 4.674182e-02 5.130127e-02 5.684684e-02 6.359080e-02 7.338302e-02 8.815700e-02 
#98%          99%         100% 
#1.123642e-01 1.781564e-01 4.423955e+00  
Osl1_cold_co <- rownames(Osl1_cold[rowVars(Osl1_cold)<0.07,])

Osl1_freeze <- expr.vst[,7:9]
head(Osl1_freeze)
quantile (rowVars(Osl1_freeze),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 9.735567e-05 
#7%           8%           9%          10%          11%          12%          13% 
#2.134899e-04 3.435740e-04 4.699275e-04 5.784735e-04 6.953117e-04 8.282822e-04 9.627201e-04 
#14%          15%          16%          17%          18%          19%          20% 
#1.082948e-03 1.213263e-03 1.349384e-03 1.487207e-03 1.633637e-03 1.771180e-03 1.894851e-03 
#21%          22%          23%          24%          25%          26%          27% 
#2.043853e-03 2.186521e-03 2.342623e-03 2.495258e-03 2.644824e-03 2.798796e-03 2.961951e-03 
#28%          29%          30%          31%          32%          33%          34% 
#3.109313e-03 3.271912e-03 3.433689e-03 3.593876e-03 3.768578e-03 3.924817e-03 4.114887e-03 
#35%          36%          37%          38%          39%          40%          41% 
#4.318983e-03 4.537250e-03 4.754058e-03 4.966156e-03 5.091101e-03 5.242030e-03 5.498596e-03 
#42%          43%          44%          45%          46%          47%          48% 
#5.715774e-03 5.972330e-03 6.199768e-03 6.446809e-03 6.720914e-03 6.964252e-03 7.240329e-03 
#49%          50%          51%          52%          53%          54%          55% 
#7.518983e-03 7.827379e-03 8.081650e-03 8.404796e-03 8.713357e-03 9.011808e-03 9.252428e-03 
#56%          57%          58%          59%          60%          61%          62% 
#9.542725e-03 9.828393e-03 1.011464e-02 1.043749e-02 1.084791e-02 1.127304e-02 1.166978e-02 
#63%          64%          65%          66%          67%          68%          69% 
#1.208866e-02 1.254895e-02 1.299276e-02 1.351627e-02 1.402521e-02 1.454275e-02 1.505944e-02 
#70%          71%          72%          73%          74%          75%          76% 
#1.558267e-02 1.617484e-02 1.677174e-02 1.751023e-02 1.819407e-02 1.912515e-02 1.991038e-02 
#77%          78%          79%          80%          81%          82%          83% 
#2.084315e-02 2.184693e-02 2.284213e-02 2.393483e-02 2.512476e-02 2.662128e-02 2.798867e-02 
#84%          85%          86%          87%          88%          89%          90% 
#2.967124e-02 3.140209e-02 3.320752e-02 3.567427e-02 3.787266e-02 4.056173e-02 4.355752e-02 
#91%          92%          93%          94%          95%          96%          97% 
#4.721651e-02 5.164927e-02 5.661962e-02 6.333995e-02 7.304742e-02 8.507682e-02 9.923588e-02 
#98%          99%         100% 
#1.247353e-01 1.846510e-01 4.514670e+00 
Osl1_freeze_co <- rownames(Osl1_freeze[rowVars(Osl1_freeze)<0.07,])

Ain1_CK <- expr.vst[,10:12]
#head(CK)
quantile (rowVars(Ain1_CK),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 7.954245e-05 1.808799e-04 2.915521e-04 3.919658e-04 5.030714e-04 
#7%           8%           9%          10%          11%          12%          13% 
#5.935945e-04 7.043022e-04 7.977162e-04 8.942801e-04 1.000654e-03 1.119705e-03 1.230420e-03 
#14%          15%          16%          17%          18%          19%          20% 
#1.345658e-03 1.475970e-03 1.583416e-03 1.697687e-03 1.820413e-03 1.939703e-03 2.066149e-03 
#21%          22%          23%          24%          25%          26%          27% 
#2.194727e-03 2.326723e-03 2.456319e-03 2.594111e-03 2.733435e-03 2.870006e-03 2.998072e-03 
#28%          29%          30%          31%          32%          33%          34% 
#3.147631e-03 3.307087e-03 3.455269e-03 3.601100e-03 3.748811e-03 3.894869e-03 4.037660e-03 
#35%          36%          37%          38%          39%          40%          41% 
#4.194382e-03 4.363861e-03 4.556065e-03 4.718872e-03 4.907671e-03 5.084400e-03 5.262033e-03 
#42%          43%          44%          45%          46%          47%          48% 
#5.450092e-03 5.656538e-03 5.858719e-03 6.050170e-03 6.235886e-03 6.447568e-03 6.670545e-03 
#49%          50%          51%          52%          53%          54%          55% 
#6.897566e-03 7.129747e-03 7.355022e-03 7.618461e-03 7.838613e-03 8.073440e-03 8.327087e-03 
#56%          57%          58%          59%          60%          61%          62% 
#8.650137e-03 8.877976e-03 9.172789e-03 9.456715e-03 9.747383e-03 1.005128e-02 1.037274e-02 
#63%          64%          65%          66%          67%          68%          69% 
#1.074022e-02 1.106129e-02 1.139955e-02 1.181215e-02 1.219201e-02 1.256299e-02 1.298600e-02 
#70%          71%          72%          73%          74%          75%          76% 
#1.343888e-02 1.387305e-02 1.437656e-02 1.485012e-02 1.532370e-02 1.586230e-02 1.646611e-02 
#77%          78%          79%          80%          81%          82%          83% 
#1.702466e-02 1.766156e-02 1.833857e-02 1.910566e-02 1.992868e-02 2.074759e-02 2.163015e-02 
#84%          85%          86%          87%          88%          89%          90% 
#2.269936e-02 2.367781e-02 2.491659e-02 2.631633e-02 2.785352e-02 2.944958e-02 3.143336e-02 
#91%          92%          93%          94%          95%          96%          97% 
#3.379516e-02 3.609615e-02 3.948931e-02 4.403239e-02 4.971209e-02 5.802245e-02 7.178744e-02 
#98%          99%         100% 
#9.485626e-02 1.526472e-01 7.015058e+00 
Ain1_CK_co <- rownames(Ain1_CK[rowVars(Ain1_CK)<0.07,]) #after check the quantile of the datasets and found 0.07 could be a good cutoff for the to 10% of the genes 
#with the highest varaince of the vst normalized counts.
Ain1_cold <- expr.vst[,13:15]
head(Ain1_cold)
quantile (rowVars(Ain1_cold),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 6.958644e-05 1.291687e-04 1.841192e-04 2.576946e-04 3.179507e-04 
#7%           8%           9%          10%          11%          12%          13% 
#3.810943e-04 4.392160e-04 5.056312e-04 5.749731e-04 6.396166e-04 7.016886e-04 7.717526e-04 
#14%          15%          16%          17%          18%          19%          20% 
#8.506895e-04 9.294631e-04 1.011149e-03 1.082002e-03 1.161141e-03 1.237996e-03 1.308495e-03 
#21%          22%          23%          24%          25%          26%          27% 
#1.392672e-03 1.472745e-03 1.564721e-03 1.641040e-03 1.725617e-03 1.823362e-03 1.912597e-03 
#28%          29%          30%          31%          32%          33%          34% 
#2.010140e-03 2.108904e-03 2.196761e-03 2.284718e-03 2.378637e-03 2.472419e-03 2.567490e-03 
#35%          36%          37%          38%          39%          40%          41% 
#2.657573e-03 2.755457e-03 2.865465e-03 2.981588e-03 3.089877e-03 3.204925e-03 3.311749e-03 
#42%          43%          44%          45%          46%          47%          48% 
#3.431216e-03 3.540298e-03 3.655345e-03 3.765416e-03 3.892610e-03 4.016357e-03 4.162309e-03 
#49%          50%          51%          52%          53%          54%          55% 
#4.287996e-03 4.417800e-03 4.548984e-03 4.694560e-03 4.824061e-03 4.948735e-03 5.092815e-03 
#56%          57%          58%          59%          60%          61%          62% 
#5.233777e-03 5.388178e-03 5.557503e-03 5.747549e-03 5.909743e-03 6.079824e-03 6.243919e-03 
#63%          64%          65%          66%          67%          68%          69% 
#6.439141e-03 6.602515e-03 6.768582e-03 6.952829e-03 7.163000e-03 7.383755e-03 7.637714e-03 
#70%          71%          72%          73%          74%          75%          76% 
#7.859409e-03 8.117770e-03 8.351428e-03 8.626116e-03 8.879072e-03 9.214568e-03 9.527571e-03 
#77%          78%          79%          80%          81%          82%          83% 
#9.817610e-03 1.012263e-02 1.049034e-02 1.093910e-02 1.136018e-02 1.177924e-02 1.222321e-02 
#84%          85%          86%          87%          88%          89%          90% 
#1.269977e-02 1.318715e-02 1.375613e-02 1.441502e-02 1.510974e-02 1.589665e-02 1.661985e-02 
#91%          92%          93%          94%          95%          96%          97% 
#1.760034e-02 1.862383e-02 2.002521e-02 2.166422e-02 2.361310e-02 2.616648e-02 3.031617e-02 
#98%          99%         100% 
#3.668188e-02 5.097373e-02 4.695292e+00
Ain1_cold_co <- rownames(Ain1_cold[rowVars(Ain1_cold)<0.07,])

Ain1_freeze <- expr.vst[,16:18]
head(Ain1_freeze)
quantile (rowVars(Ain1_freeze),probs=seq(0,1,0.01))
#0%           1%           2%           3%           4%           5%           6% 
#0.000000e+00 0.000000e+00 8.006294e-05 1.867370e-04 2.901763e-04 3.956246e-04 5.128187e-04 
#7%           8%           9%          10%          11%          12%          13% 
#6.200345e-04 7.370504e-04 8.512390e-04 9.728242e-04 1.084524e-03 1.200879e-03 1.314592e-03 
#14%          15%          16%          17%          18%          19%          20% 
#1.437402e-03 1.563054e-03 1.687488e-03 1.817352e-03 1.933529e-03 2.057329e-03 2.203186e-03 
#21%          22%          23%          24%          25%          26%          27% 
#2.334130e-03 2.478747e-03 2.623274e-03 2.775526e-03 2.936783e-03 3.093676e-03 3.232546e-03 
#28%          29%          30%          31%          32%          33%          34% 
#3.378194e-03 3.529306e-03 3.663283e-03 3.821470e-03 3.970507e-03 4.126395e-03 4.291864e-03 
#35%          36%          37%          38%          39%          40%          41% 
#4.461648e-03 4.643328e-03 4.804807e-03 5.006942e-03 5.215817e-03 5.379292e-03 5.572892e-03 
#42%          43%          44%          45%          46%          47%          48% 
#5.791041e-03 5.979049e-03 6.193531e-03 6.415777e-03 6.621985e-03 6.819792e-03 7.009321e-03 
#49%          50%          51%          52%          53%          54%          55% 
#7.232256e-03 7.452206e-03 7.684479e-03 7.902139e-03 8.141342e-03 8.399668e-03 8.652774e-03 
#56%          57%          58%          59%          60%          61%          62% 
#8.885620e-03 9.154633e-03 9.368305e-03 9.691239e-03 1.002505e-02 1.035521e-02 1.066266e-02 
#63%          64%          65%          66%          67%          68%          69% 
#1.100807e-02 1.132110e-02 1.166088e-02 1.200579e-02 1.243222e-02 1.282162e-02 1.318798e-02 
#70%          71%          72%          73%          74%          75%          76% 
#1.368704e-02 1.412728e-02 1.462285e-02 1.506723e-02 1.559384e-02 1.611842e-02 1.668053e-02 
#77%          78%          79%          80%          81%          82%          83% 
#1.728931e-02 1.795392e-02 1.862774e-02 1.936583e-02 2.011023e-02 2.101415e-02 2.197787e-02 
#84%          85%          86%          87%          88%          89%          90% 
#2.299664e-02 2.423397e-02 2.539145e-02 2.665646e-02 2.825018e-02 2.972921e-02 3.175234e-02 
#91%          92%          93%          94%          95%          96%          97% 
#3.416447e-02 3.704713e-02 4.041881e-02 4.422444e-02 4.984432e-02 5.732636e-02 6.929996e-02 
#98%          99%         100% 
#8.844654e-02 1.363025e-01 6.733631e+00 
Ain1_freeze_co <- rownames(Ain1_freeze[rowVars(Ain1_freeze)<0.07,])

#all of treatment has variance <0.07
low_var_gene <- intersect(intersect(intersect(intersect(intersect(Osl1_CK_co,Osl1_cold_co),Osl1_freeze_co),Ain1_CK_co),Ain1_cold_co),Ain1_freeze_co)# here are the genes with loc varaince for each treatment 

length(low_var_gene)
#19691
counts <- counts[rownames(counts) %in% low_var_gene,]
nrow(counts)
#19691
grep("Brasy1G039900", rownames(counts))#check if this CBF got filtered
#226
#it suggested that the filtering seems good!
#19691/36927=0.5332413 only keep around 60% genes

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
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/normalized_VST_final_cold_formal#.txt",
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
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/PCA_cold_syl_novargens.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
plotPCA(rld, intgroup="group")
dev.off()

head(rld)
library( "genefilter" )
#install.packages("gplots")
library(gplots)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 1155 )
head(topVarGenes)
#dev.off()
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_1155_DGE_cold_syl_novargens.pdf",width = 15,height = 15)

par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row', 
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", cold="#11C5F6", freeze="#1149F6")[
             colData(rld)$treatment ] )
dev.off()
resultsNames(ddsTC)
##control versus treatments for Osl1
res_Ain1_cold <- results(ddsTC, contrast=list("groupAin1.cold", "groupAin1.CK"))
EnhancedVolcano(res_Osl1_cold,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Ain1_cold',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=2,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/cold_syl_ain1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#colnames(ddsTC)
res_Ain1_freeze <- results(ddsTC, contrast=list("groupAin1.freeze", "groupAin1.CK"))
EnhancedVolcano(res_Ain1_freeze,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Ain1_freeze',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=2,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/freeze_syl_ain1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)

res_Osl1_cold <- results(ddsTC, contrast=list("groupOsl1.cold", "groupOsl1.CK"))
EnhancedVolcano(res_Osl1_cold,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Osl1_cold',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=2,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/cold_syl_Osl1.pdf",
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 22,
  units = "cm",
  dpi = 300,
)
#colnames(ddsTC)
res_Osl1_freeze <- results(ddsTC, contrast=list("groupOsl1.freeze", "groupOsl1.CK"))
EnhancedVolcano(res_Osl1_freeze,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                title =  'Osl1_freeze',
                xlim = c(-25, 25),
                ylim = c(0,320),
                FCcutoff=2,
                pCutoff = 0.05,
                pointSize = 1.5,
                colAlpha = 0.5)
ggsave(
  "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/freeze_syl_Osl1.pdf",
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
res_Ain1_cold <- res_Ain1_cold[!is.na(res_Ain1_cold$padj),]
sigres_res_Ain1_cold_up <- res_Ain1_cold[(res_Ain1_cold$padj<0.05 & res_Ain1_cold$log2FoldChange>= 1),]
nb_Ain1_cold_up <- nrow(sigres_res_Ain1_cold_up)
#1057 delte the big variant genes
write.table(x = sigres_res_Ain1_cold_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Ain1_cold_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Ain1_cold_down <- res_Ain1_cold[(res_Ain1_cold$padj<0.05 & res_Ain1_cold$log2FoldChange<= -1),]
nb_Ain1_cold_down <- nrow(sigres_res_Ain1_cold_down)
#1623
write.table(x = sigres_res_Ain1_cold_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Ain1_cold_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Ain1_cold_DEG <- res_Ain1_cold[(res_Ain1_cold$padj<0.05 & abs(res_Ain1_cold$log2FoldChange) >= 1 ),]
nb_sigres_res_Ain1_cold_DEG <- nrow(sigres_res_Ain1_cold_DEG)
#2680
#set the logfoldchanges as 5 and extract those genes
sigres_res_Ain1_cold_DEG_strict <- res_Ain1_cold[(res_Ain1_cold$padj<0.05 & abs(res_Ain1_cold$log2FoldChange) >= 2),]
head(sigres_res_Ain1_cold_DEG_strict)
nrow(sigres_res_Ain1_cold_DEG_strict)
#507
write.table(x = sigres_res_Ain1_cold_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Ain1_cold_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##freeze
res_Ain1_freeze <- res_Ain1_freeze[!is.na(res_Ain1_freeze$padj),]
sigres_res_Ain1_freeze_up <- res_Ain1_freeze[(res_Ain1_freeze$padj<0.05 & res_Ain1_freeze$log2FoldChange>= 1),]
nb_Ain1_freeze_up <- nrow(sigres_res_Ain1_freeze_up)
#1456
write.table(x = sigres_res_Ain1_freeze_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Ain1_freeze_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Ain1_freeze_down <- res_Ain1_freeze[(res_Ain1_freeze$padj<0.05 & res_Ain1_freeze$log2FoldChange<= -1),]
nb_Ain1_freeze_down <- nrow(sigres_res_Ain1_freeze_down)
#1831
write.table(x = sigres_res_Ain1_freeze_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Ain1_freeze_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Ain1_freeze_DEG <- res_Ain1_freeze[(res_Ain1_freeze$padj<0.05 & abs(res_Ain1_freeze$log2FoldChange) >= 1 ),]
nb_sigres_res_Ain1_freeze_DEG <- nrow(sigres_res_Ain1_freeze_DEG)
#3287
#set the logfoldchanges as 2 and extract those genes
sigres_res_Ain1_freeze_DEG_strict <- res_Ain1_freeze[(res_Ain1_freeze$padj<0.05 & abs(res_Ain1_freeze$log2FoldChange) >= 2),]
head(sigres_res_Ain1_freeze_DEG_strict)
nrow(sigres_res_Ain1_freeze_DEG_strict)
#721
write.table(x = sigres_res_Ain1_freeze_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Ain1_freeze_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

#cold
res_Osl1_cold <- res_Osl1_cold[!is.na(res_Osl1_cold$padj),]
sigres_res_Osl1_cold_up <- res_Osl1_cold[(res_Osl1_cold$padj<0.05 & res_Osl1_cold$log2FoldChange>= 1),]
nb_Osl1_cold_up <- nrow(sigres_res_Osl1_cold_up)
#1356 delte the big variant genes
write.table(x = sigres_res_Osl1_cold_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Osl1_cold_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Osl1_cold_down <- res_Osl1_cold[(res_Osl1_cold$padj<0.05 & res_Osl1_cold$log2FoldChange<= -1),]
nb_Osl1_cold_down <- nrow(sigres_res_Osl1_cold_down)
#1700
write.table(x = sigres_res_Osl1_cold_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Osl1_cold_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Osl1_cold_DEG <- res_Osl1_cold[(res_Osl1_cold$padj<0.05 & abs(res_Osl1_cold$log2FoldChange) >= 1 ),]
nb_sigres_res_Osl1_cold_DEG <- nrow(sigres_res_Osl1_cold_DEG)
#3056
#set the logfoldchanges as 5 and extract those genes
sigres_res_Osl1_cold_DEG_strict <- res_Osl1_cold[(res_Osl1_cold$padj<0.05 & abs(res_Osl1_cold$log2FoldChange) >= 2),]
head(sigres_res_Osl1_cold_DEG_strict)
nrow(sigres_res_Osl1_cold_DEG_strict)
#589
write.table(x = sigres_res_Osl1_cold_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Osl1_cold_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)
##freeze
res_Osl1_freeze <- res_Osl1_freeze[!is.na(res_Osl1_freeze$padj),]
sigres_res_Osl1_freeze_up <- res_Osl1_freeze[(res_Osl1_freeze$padj<0.05 & res_Osl1_freeze$log2FoldChange>= 1),]
nb_Osl1_freeze_up <- nrow(sigres_res_Osl1_freeze_up)
#11891
write.table(x = sigres_res_Osl1_freeze_up,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Osl1_freeze_up.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Osl1_freeze_down <- res_Osl1_freeze[(res_Osl1_freeze$padj<0.05 & res_Osl1_freeze$log2FoldChange<= -1),]
nb_Osl1_freeze_down <- nrow(sigres_res_Osl1_freeze_down)
#2624
write.table(x = sigres_res_Osl1_freeze_down,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Osl1_freeze_down.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)

sigres_res_Osl1_freeze_DEG <- res_Osl1_freeze[(res_Osl1_freeze$padj<0.05 & abs(res_Osl1_freeze$log2FoldChange) >= 1 ),]
nb_sigres_res_Osl1_freeze_DEG <- nrow(sigres_res_Osl1_freeze_DEG)
#4515
#set the logfoldchanges as 2 and extract those genes
sigres_res_Osl1_freeze_DEG_strict <- res_Osl1_freeze[(res_Osl1_freeze$padj<0.05 & abs(res_Osl1_freeze$log2FoldChange) >= 2),]
head(sigres_res_Osl1_freeze_DEG_strict)
nrow(sigres_res_Osl1_freeze_DEG_strict)
#1293
write.table(x = sigres_res_Osl1_freeze_DEG,
            file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/sigres_res_Osl1_freeze_DEG.txt",
            sep = "\t",
            eol = "\n",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE
)


#top gene list:
top_gene_list_str <- union(union(union(rownames(sigres_res_Ain1_freeze_DEG_strict),rownames(sigres_res_Ain1_cold_DEG_strict)),rownames(sigres_res_Osl1_cold_DEG_strict)),rownames(sigres_res_Osl1_freeze_DEG_strict))
length(top_gene_list_str)
#1959
top_gene_list <- union(union(union(rownames(sigres_res_Ain1_freeze_DEG),rownames(sigres_res_Ain1_cold_DEG)),rownames(sigres_res_Osl1_cold_DEG)),rownames(sigres_res_Osl1_freeze_DEG))
length(top_gene_list)
#7136

#plot the genes with logfolderchanges >1
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_5155_DGE_cold_syl_ain_osl.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[row.names(assay(rld)) %in% top_gene_list,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", cold="#11C5F6", freeze="#1149F6")[
             colData(rld)$treatment ] )
dev.off()
#plot the genes with logfolderchanges >2
pdf(file = "/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/heatmap_1959_DGE_cold_syl_ain_osl.pdf",width = 15,height = 15)
par(mar=c(1,1,1,1))
#heatmap(assay(rld)[ topVarGenes, ],
#        labRow = "",
#        scale = "row")
palette <- colorRampPalette(c("red","white","blue"))(256)

heatmap.2( assay(rld)[row.names(assay(rld)) %in% top_gene_list_str,], scale="row", Colv=FALSE, 
           trace="none", dendrogram='row',
           #Rowv=FALSE,
           #Colv=FALSE,
           col = palette,
           lhei = c(0.8,7),
           ColSideColors = c( CK="gray", cold="#11C5F6", freeze="#1149F6")[
             colData(rld)$treatment ] )
dev.off()

#
library(ggplot2)
library(plyr)
library(reshape)
library(scales)
library(tidyverse)
library(forcats)
cold_nb <- data.frame(Accession=rep(c("Ain1","Osl1"),each=4),
                      Categories=rep(c("Up", "Down"), each=2),
                      Treatments=rep(c("cold", "freeze"),2),
                      Number=c(nb_Ain1_cold_up,nb_Ain1_freeze_up,nb_Ain1_cold_down,nb_Ain1_freeze_down,nb_Osl1_cold_up,nb_Osl1_freeze_up,nb_Osl1_cold_down,nb_Osl1_freeze_down))

ggplot(cold_nb, aes(x=factor(Treatments),fill=Accession )) + 
  geom_bar(data = subset(cold_nb, Categories == "Up"),
           aes(y=Number, fill = Accession),stat = "identity", position ="dodge")+
  scale_fill_manual(values = c("black","grey"))+
  theme(axis.ticks.length=unit(.25, "cm"),
        axis.text.y = element_text(color = "black", size = 20, angle = 0, hjust = 1, vjust = 0, face = "bold"),
        axis.text.x = element_text(color = "black", size = 25, angle = 0, hjust = 1, vjust = 0, face = "bold"))+
  ylim(-3000, 3000)+ theme_bw()

last_plot() + geom_bar(data = subset(cold_nb, Categories == "Down"), 
                       aes(y= -Number, fill = Accession), stat = 'identity', position = "dodge") +
  xlab("") + 
  ylab("Down  -  Up")+
  geom_hline(yintercept = 0,colour = "grey90") 

ggsave("/global/projectb/scratch/llei2019/RNAseq_Syl_sgordon/cold_DEG_nb_accessions.pdf", width = 15, height = 12)

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
