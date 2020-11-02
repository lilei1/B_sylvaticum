#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
library(tidyr)
#########################################################
### B) Reading in data and transform it into matrix format
#########################################################

#data <- read.csv("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/unique_heat_down_dict_GO_1h_2h_5h_10h_24h.txt",row.names = 3,sep="\t",header = T)
data <- read.csv("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/unique_drought_up_dict_GO_1h_2h_5h_10h_24h.txt",row.names = 3,sep="\t",header = T)
#data <- read.csv("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/unique_drought_down_dict_GO_1h_2h_5h_10h_24h.txt",row.names = 3,sep="\t",header = T)
#data <- read.csv("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/unique_salt_up_dict_GO_1h_2h_5h_10h_24h.txt",row.names = 3,sep="\t",header = T)
#data <- read.csv("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/unique_salt_down_dict_GO_2h_5h_10h_24h.txt",row.names = 3,sep="\t",header = T)

head(data)
subdata <- data[,c("drought1h","drought2h","drought5h","drought10h","drought24h")] 
#subdata <- data[,c("salt1h","salt2h","salt5h","salt10h","salt24h")] 
mat_data <- data.matrix(subdata)  # transform column 2-5 into a matrix
single_column <- gather(subdata,key = GO, value = foldchanges,na.rm =T)
head(single_column)
quantile(single_column$foldchanges,probs = seq(0,1,0.01))
#down_salt
#0%      1%      2%      3%      4%      5%      6%      7%      8%      9%     10%     11%     12%     13%     14%     15%     16% 
#1.1500  1.2022  1.2800  1.4250  1.4900  1.4900  1.5092  1.5312  1.5428  1.5632  1.5980  1.6252  1.6484  1.6554  1.6660  1.6950  1.7352 
#17%     18%     19%     20%     21%     22%     23%     24%     25%     26%     27%     28%     29%     30%     31%     32%     33% 
#1.7874  1.8044  1.8134  1.9120  1.9818  1.9876  2.0682  2.2016  2.2250  2.2316  2.2432  2.2548  2.2664  2.2820  2.2994  2.3168  2.3454 
#34%     35%     36%     37%     38%     39%     40%     41%     42%     43%     44%     45%     46%     47%     48%     49%     50% 
#2.4092  2.4430  2.4488  2.5144  2.5916  2.6148  2.6300  2.6300  2.6552  2.6958  2.7624  2.8230  2.8404  2.8578  2.8752  2.8842  2.8900 
#51%     52%     53%     54%     55%     56%     57%     58%     59%     60%     61%     62%     63%     64%     65%     66%     67% 
#2.9306  2.9664  2.9896  3.0480  3.1350  3.2604  3.3818  3.3992  3.4430  3.5300  3.5600  3.5600  3.6140  3.6684  3.7090  3.7300  3.7300 
#68%     69%     70%     71%     72%     73%     74%     75%     76%     77%     78%     79%     80%     81%     82%     83%     84% 
#3.7300  3.7334  3.8320  3.9468  4.0976  4.1600  4.1600  4.2150  4.2852  4.3954  4.4768  4.5174  4.5740  4.6378  4.8416  5.0476  5.2448 
#85%     86%     87%     88%     89%     90%     91%     92%     93%     94%     95%     96%     97%     98%     99%    100% 
#  5.3610  5.4016  5.4100  5.4368  5.8254  6.0800  6.0800  6.5120  7.2080  8.6944 10.0910 10.6188 11.1050 11.5400 12.7352 14.2200 

#up_salt
#0%      1%      2%      3%      4%      5%      6%      7%      8%      9%     10%     11%     12%     13%     14%     15%     16% 
#1.1400  1.1500  1.1516  1.1600  1.1632  1.1740  1.1848  1.2012  1.2164  1.2416  1.2580  1.2600  1.2792  1.2808  1.3012  1.3140  1.3328 
#17%     18%     19%     20%     21%     22%     23%     24%     25%     26%     27%     28%     29%     30%     31%     32%     33% 
#1.3436  1.3544  1.3600  1.3660  1.3904  1.4000  1.4000  1.4000  1.4100  1.4108  1.4200  1.4248  1.4400  1.4440  1.4500  1.4612  1.4700 
#34%     35%     36%     37%     38%     39%     40%     41%     42%     43%     44%     45%     46%     47%     48%     49%     50% 
#1.4772  1.4880  1.4988  1.5000  1.5000  1.5000  1.5040  1.5256  1.5400  1.5400  1.5400  1.5640  1.5868  1.5900  1.5900  1.5992  1.6100 
#51%     52%     53%     54%     55%     56%     57%     58%     59%     60%     61%     62%     63%     64%     65%     66%     67% 
#1.6124  1.6432  1.6672  1.6900  1.6900  1.6948  1.7000  1.7000  1.7144  1.7600  1.7788  1.7896  1.8020  1.8524  1.8760  1.9028  1.9100 
#68%     69%     70%     71%     72%     73%     74%     75%     76%     77%     78%     79%     80%     81%     82%     83%     84% 
#1.9628  2.0352  2.0940  2.1436  2.1576  2.1600  2.2244  2.3700  2.3900  2.4172  2.5744  2.6392  2.6800  2.9968  3.4800  3.9036  4.0872 
#85%     86%     87%     88%     89%     90%     91%     92%     93%     94%     95%     96%     97%     98%     99%    100% 
#  4.1140  4.2520  4.2700  4.2712  4.3588  4.7920  4.8000  4.8000  4.8000  4.8312  4.8840  4.9000  4.9532  5.3396  5.4100 10.7000 

#down-drought
#0%    1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13%   14%   15%   16%   17%   18%   19%   20%   21%   22% 
#1.070 1.090 1.120 1.130 1.190 1.205 1.220 1.225 1.230 1.230 1.230 1.240 1.250 1.250 1.260 1.260 1.270 1.280 1.280 1.295 1.300 1.300 1.300 
#23%   24%   25%   26%   27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39%   40%   41%   42%   43%   44%   45% 
#1.310 1.320 1.325 1.330 1.330 1.330 1.330 1.330 1.335 1.340 1.340 1.350 1.365 1.370 1.370 1.370 1.370 1.380 1.380 1.380 1.390 1.390 1.395 
#46%   47%   48%   49%   50%   51%   52%   53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65%   66%   67%   68% 
#1.400 1.405 1.410 1.420 1.420 1.425 1.460 1.480 1.490 1.500 1.500 1.575 1.600 1.685 1.730 1.735 1.740 1.800 1.830 1.880 1.950 1.980 2.010 
#69%   70%   71%   72%   73%   74%   75%   76%   77%   78%   79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
#2.035 2.070 2.110 2.160 2.205 2.320 2.390 2.420 2.470 2.470 2.470 2.510 2.545 2.640 2.710 3.010 3.160 3.280 3.380 3.410 3.555 3.590 3.825 
#92%   93%   94%   95%   96%   97%   98%   99%  100% 
#3.910 4.325 4.550 4.680 4.740 5.315 5.530 6.245 6.630 
#drought-up
#0%    1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13%   14%   15%   16%   17%   18%   19%   20%   21%   22% 
#1.070 1.110 1.120 1.150 1.150 1.160 1.202 1.219 1.236 1.263 1.270 1.284 1.294 1.311 1.320 1.335 1.340 1.340 1.356 1.363 1.370 1.380 1.380 
#23%   24%   25%   26%   27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39%   40%   41%   42%   43%   44%   45% 
#1.380 1.388 1.410 1.410 1.410 1.416 1.430 1.430 1.430 1.440 1.450 1.450 1.465 1.470 1.470 1.516 1.530 1.530 1.547 1.568 1.581 1.590 1.605 
#46%   47%   48%   49%   50%   51%   52%   53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65%   66%   67%   68% 
#1.610 1.610 1.620 1.620 1.620 1.637 1.644 1.660 1.660 1.660 1.664 1.689 1.702 1.740 1.740 1.740 1.852 1.882 1.916 1.945 1.976 2.018 2.070 
#69%   70%   71%   72%   73%   74%   75%   76%   77%   78%   79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
#2.116 2.130 2.168 2.200 2.205 2.274 2.305 2.362 2.410 2.410 2.410 2.480 2.525 2.588 2.601 2.618 2.625 2.822 2.848 2.860 2.880 2.880 3.056 
#92%   93%   94%   95%   96%   97%   98%   99%  100% 
#3.180 3.208 3.280 3.280 3.388 3.639 4.050 6.764 7.520 
#> 
# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)
my_palette <- colorRampPalette(c("white","yellow", "red"))(n = 299)
col_breaks = c(seq(0,1.45,length=100),  # for red
               seq(1.46,1.97,length=100),           # for yellow
               seq(1.98,7.52,length=100)) 

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.47,length=100),  # for red
#               seq(1.48,1.90,length=100),           # for yellow
#               seq(1.91,10.7,length=100))             # for green
#col_breaks = c(seq(0,2.34,length=100),  # for red
 #              seq(2.35,3.73,length=100),           # for yellow
 #              seq(3.74,14.22,length=100))             # for green
 #creates a 5 x 5 inch image
pdf("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/heatmaps_drought_up.pdf",    # create PNG for the heat map        
    width = 15,        # 5 x 300 pixels
    height =20
)        # smaller font size
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Foldchange", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",# Don't draw dendrogram
          distfun=dist_no_na,
          cexRow=2,
          cexCol=2,
          Colv="NA",
          key =T,
          lwid = c(0.6, 4, 0.6),
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ),
          lhei = c(0.43, 2.6, 0.6))            # turn off column clustering

dev.off()               # close the PNG device



#Cellular
cellular <- data[(data$cat == "C"),]
#C_subdata <- cellular[,c("heat1h","heat2h","heat5h","heat10h","heat24h")] 
#C_subdata <- cellular[,c("salt1h","salt2h","salt5h","salt10h","salt24h")] 
C_subdata <- cellular[,c("drought1h","drought2h","drought5h","drought10h","drought24h")] 

#subdata <- data[,c("drought1h","drought2h","drought5h","drought10h","drought24h")] 

head(C_subdata)
#head(subdata)
# assign labels in column 1 to "rnames"
mat_data <- data.matrix(C_subdata)  # transform column 2-5 into a matrix
#mat_data <- data.matrix(subdata)  # transform column 2-5 into a matrix
#rownames(mat_data) <- rnames 
head(mat_data)
#single_column <- gather(C_subdata,key = GO, value = foldchanges,na.rm =T)
#single_column <- gather(subdata,key = GO, value = foldchanges,na.rm =T)

#head(single_column)
#quantile(single_column$foldchanges,probs = seq(0,1,0.01))
#0%        1%        2%        3%        4%        5%        6%        7%        8%        9%       10% 
#1.303440  1.304030  1.304620  1.305210  1.320673  1.338261  1.355848  1.360245  1.360245  1.360245  1.372196 
#11%       12%       13%       14%       15%       16%       17%       18%       19%       20%       21% 
#1.391318  1.410439  1.420000  1.420000  1.420000  1.448800  1.525600  1.602400  1.660000  1.660000  1.660000 
#22%       23%       24%       25%       26%       27%       28%       29%       30%       31%       32% 
#1.668800  1.739200  1.809600  1.880000  1.880000  1.880000  1.880000  1.891200  1.904000  1.916800  1.970400 
#33%       34%       35%       36%       37%       38%       39%       40%       41%       42%       43% 
#2.037600  2.104800  2.452000  2.967200  3.482400  3.794400  3.903200  4.012000  4.088007  4.109360  4.130713 
#44%       45%       46%       47%       48%       49%       50%       51%       52%       53%       54% 
#4.173389  4.280036  4.386684  4.480800  4.487200  4.493600  4.500000  4.615200  4.730400  4.845600  4.882400 
#55%       56%       57%       58%       59%       60%       61%       62%       63%       64%       65% 
#4.908000  4.933600  5.026400  5.141600  5.256800  5.516000  5.861600  6.207200  6.439200  6.557600  6.676000 
#66%       67%       68%       69%       70%       71%       72%       73%       74%       75%       76% 
#6.753600  6.763200  6.772800  6.807251  6.916257  7.025263  7.141817  7.311212  7.480606  7.650000  7.650000 
#77%       78%       79%       80%       81%       82%       83%       84%       85%       86%       87% 
#7.650000  7.650000  7.669600  7.692000  7.714400  7.787200  7.876800  7.966400  8.236000  8.613600  8.991200 
#88%       89%       90%       91%       92%       93%       94%       95%       96%       97%       98% 
#9.311200  9.573600  9.836000 10.025200 10.092400 10.159600 10.413200 11.226000 12.038800 12.969200 14.722800 
#99%      100% 
#  16.476400 18.230000 

#down
#0%      1%      2%      3%      4%      5%      6%      7%      8%      9%     10%     11%     12%     13% 
#  1.1000  1.1500  1.1600  1.1651  1.1768  1.1800  1.1900  1.2119  1.2336  1.2600  1.2600  1.2700  1.2800  1.2900 
#14%     15%     16%     17%     18%     19%     20%     21%     22%     23%     24%     25%     26%     27% 
#  1.2900  1.3000  1.3000  1.3000  1.3106  1.3200  1.3240  1.3500  1.3500  1.3782  1.3908  1.4000  1.4100  1.4159 
#28%     29%     30%     31%     32%     33%     34%     35%     36%     37%     38%     39%     40%     41% 
#  1.4200  1.4293  1.4400  1.4500  1.4644  1.4861  1.4978  1.5100  1.5200  1.5358  1.5700  1.5863  1.6100  1.6100 
#42%     43%     44%     45%     46%     47%     48%     49%     50%     51%     52%     53%     54%     55% 
#  1.6314  1.6431  1.6548  1.6700  1.6882  1.7099  1.7300  1.7400  1.7500  1.7667  1.7884  1.8802  1.9126  1.9870 
#56%     57%     58%     59%     60%     61%     62%     63%     64%     65%     66%     67%     68%     69% 
#  2.0152  2.0369  2.1786  2.1900  2.2300  2.2948  2.3254  2.3500  2.3964  2.4305  2.4754  2.5300  2.5412  2.5846 
#70%     71%     72%     73%     74%     75%     76%     77%     78%     79%     80%     81%     82%     83% 
#  2.6200  2.6514  2.7324  2.8441  2.8616  2.9450  3.0952  3.1627  3.2000  3.2644  3.3800  3.4000  3.5552  3.7411 
#84%     85%     86%     87%     88%     89%     90%     91%     92%     93%     94%     95%     96%     97% 
#  3.7824  3.8400  4.0386  4.1411  4.4204  4.6069  4.9430  5.0300  5.0300  5.1200  5.1700  5.6170  6.8796  7.2400 
#98%     99%    100% 
#7.5962 10.2327 12.0700 

#drought-up
#0%    1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13%   14%   15%   16%   17%   18%   19%   20%   21%   22% 
#1.070 1.110 1.120 1.150 1.150 1.160 1.202 1.219 1.236 1.263 1.270 1.284 1.294 1.311 1.320 1.335 1.340 1.340 1.356 1.363 1.370 1.380 1.380 
#23%   24%   25%   26%   27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39%   40%   41%   42%   43%   44%   45% 
#1.380 1.388 1.410 1.410 1.410 1.416 1.430 1.430 1.430 1.440 1.450 1.450 1.465 1.470 1.470 1.516 1.530 1.530 1.547 1.568 1.581 1.590 1.605 
#46%   47%   48%   49%   50%   51%   52%   53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65%   66%   67%   68% 
#1.610 1.610 1.620 1.620 1.620 1.637 1.644 1.660 1.660 1.660 1.664 1.689 1.702 1.740 1.740 1.740 1.852 1.882 1.916 1.945 1.976 2.018 2.070 
#69%   70%   71%   72%   73%   74%   75%   76%   77%   78%   79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
#2.116 2.130 2.168 2.200 2.205 2.274 2.305 2.362 2.410 2.410 2.410 2.480 2.525 2.588 2.601 2.618 2.625 2.822 2.848 2.860 2.880 2.880 3.056 
#92%   93%   94%   95%   96%   97%   98%   99%  100% 
#3.180 3.208 3.280 3.280 3.388 3.639 4.050 6.764 7.520 

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
my_palette <- colorRampPalette(c("white","yellow", "red"))(n = 299)
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)
col_breaks = c(seq(0,1.45,length=100),  # for red
               seq(1.46,1.97,length=100),           # for yellow
               seq(1.98,7.52,length=100)) 
# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.47,length=100),  # for red
#               seq(1.48,1.90,length=100),           # for yellow
#               seq(1.91,10.7,length=100))             # for green
#col_breaks = c(seq(0,2.34,length=100),  # for red
#               seq(2.35,3.73,length=100),           # for yellow
#               seq(3.74,14.22,length=100)) 
# creates a 5 x 5 inch image
pdf("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/heatmaps_drought_up_C.pdf",    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height =10
    )        # smaller font size
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Foldchange", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",# Don't draw dendrogram
          distfun=dist_no_na,
          cexRow=2,
          cexCol=2,
          Colv="NA",
          key =T,
          lwid = c(0.6, 4, 0.6),
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ),
          lhei = c(0.43, 2.6, 0.6))            # turn off column clustering

dev.off()               # close the PNG device

#F
F <- data[(data$cat == "F"),]
#F_subdata <- F[,c("salt1h","salt2h","salt5h","salt10h","salt24h")] 
F_subdata <- F[,c("drought1h","drought2h","drought5h","drought10h","drought24h")] 

head(F_subdata)
# assign labels in column 1 to "rnames"
mat_data <- data.matrix(F_subdata)  # transform column 2-5 into a matrix
#rownames(mat_data) <- rnames 
head(mat_data)
single_column <- gather(subdata,key = GO, value = foldchanges,na.rm =T)
#head(single_column)
#quantile(single_column$foldchanges,probs = seq(0,1,0.01))
#0%        1%        2%        3%        4%        5%        6%        7%        8%        9%       10% 
#1.303440  1.304030  1.304620  1.305210  1.320673  1.338261  1.355848  1.360245  1.360245  1.360245  1.372196 
#11%       12%       13%       14%       15%       16%       17%       18%       19%       20%       21% 
#1.391318  1.410439  1.420000  1.420000  1.420000  1.448800  1.525600  1.602400  1.660000  1.660000  1.660000 
#22%       23%       24%       25%       26%       27%       28%       29%       30%       31%       32% 
#1.668800  1.739200  1.809600  1.880000  1.880000  1.880000  1.880000  1.891200  1.904000  1.916800  1.970400 
#33%       34%       35%       36%       37%       38%       39%       40%       41%       42%       43% 
#2.037600  2.104800  2.452000  2.967200  3.482400  3.794400  3.903200  4.012000  4.088007  4.109360  4.130713 
#44%       45%       46%       47%       48%       49%       50%       51%       52%       53%       54% 
#4.173389  4.280036  4.386684  4.480800  4.487200  4.493600  4.500000  4.615200  4.730400  4.845600  4.882400 
#55%       56%       57%       58%       59%       60%       61%       62%       63%       64%       65% 
#4.908000  4.933600  5.026400  5.141600  5.256800  5.516000  5.861600  6.207200  6.439200  6.557600  6.676000 
#66%       67%       68%       69%       70%       71%       72%       73%       74%       75%       76% 
#6.753600  6.763200  6.772800  6.807251  6.916257  7.025263  7.141817  7.311212  7.480606  7.650000  7.650000 
#77%       78%       79%       80%       81%       82%       83%       84%       85%       86%       87% 
#7.650000  7.650000  7.669600  7.692000  7.714400  7.787200  7.876800  7.966400  8.236000  8.613600  8.991200 
#88%       89%       90%       91%       92%       93%       94%       95%       96%       97%       98% 
#9.311200  9.573600  9.836000 10.025200 10.092400 10.159600 10.413200 11.226000 12.038800 12.969200 14.722800 
#99%      100% 
#  16.476400 18.230000 

#down
#0%      1%      2%      3%      4%      5%      6%      7%      8%      9%     10%     11%     12%     13% 
#  1.1000  1.1500  1.1600  1.1651  1.1768  1.1800  1.1900  1.2119  1.2336  1.2600  1.2600  1.2700  1.2800  1.2900 
#14%     15%     16%     17%     18%     19%     20%     21%     22%     23%     24%     25%     26%     27% 
#  1.2900  1.3000  1.3000  1.3000  1.3106  1.3200  1.3240  1.3500  1.3500  1.3782  1.3908  1.4000  1.4100  1.4159 
#28%     29%     30%     31%     32%     33%     34%     35%     36%     37%     38%     39%     40%     41% 
#  1.4200  1.4293  1.4400  1.4500  1.4644  1.4861  1.4978  1.5100  1.5200  1.5358  1.5700  1.5863  1.6100  1.6100 
#42%     43%     44%     45%     46%     47%     48%     49%     50%     51%     52%     53%     54%     55% 
#  1.6314  1.6431  1.6548  1.6700  1.6882  1.7099  1.7300  1.7400  1.7500  1.7667  1.7884  1.8802  1.9126  1.9870 
#56%     57%     58%     59%     60%     61%     62%     63%     64%     65%     66%     67%     68%     69% 
#  2.0152  2.0369  2.1786  2.1900  2.2300  2.2948  2.3254  2.3500  2.3964  2.4305  2.4754  2.5300  2.5412  2.5846 
#70%     71%     72%     73%     74%     75%     76%     77%     78%     79%     80%     81%     82%     83% 
#  2.6200  2.6514  2.7324  2.8441  2.8616  2.9450  3.0952  3.1627  3.2000  3.2644  3.3800  3.4000  3.5552  3.7411 
#84%     85%     86%     87%     88%     89%     90%     91%     92%     93%     94%     95%     96%     97% 
#  3.7824  3.8400  4.0386  4.1411  4.4204  4.6069  4.9430  5.0300  5.0300  5.1200  5.1700  5.6170  6.8796  7.2400 
#98%     99%    100% 
#7.5962 10.2327 12.0700 

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.48,length=100),  # for red
#               seq(1.49,2.47,length=100),           # for yellow
#               seq(2.48,12.07,length=100))             # for green
#my_palette <- colorRampPalette(c("white","yellow", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.47,length=100),  # for red
#               seq(1.48,1.90,length=100),           # for yellow
#               seq(1.91,10.7,length=100))             # for green
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)

#col_breaks = c(seq(0,2.34,length=100),  # for red
#               seq(2.35,3.73,length=100),           # for yellow
#               seq(3.74,14.22,length=100)) 
my_palette <- colorRampPalette(c("white","yellow", "red"))(n = 299)
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)
col_breaks = c(seq(0,1.45,length=100),  # for red
               seq(1.46,1.97,length=100),           # for yellow
               seq(1.98,7.52,length=100)) 

# creates a 5 x 5 inch image
pdf("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/heatmaps_drought_up_F.pdf",    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height =10
)        # smaller font size
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Foldchange", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",# Don't draw dendrogram
          distfun=dist_no_na,
          cexRow=2,
          cexCol=2,
          Colv="NA",
          key =T,
          lwid = c(0.6, 4, 0.6),
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ),
          lhei = c(0.43, 2.6, 0.6))            # turn off column clustering

dev.off()               # close the PNG device

#P
P <- data[( data$cat == "P"),]
#P_subdata <- P[,c("salt1h","salt2h","salt5h","salt10h","salt24h")] 
P_subdata <- P[,c("drought1h","drought2h","drought5h","drought10h","drought24h")] 

head(P_subdata)
# assign labels in column 1 to "rnames"
mat_data <- data.matrix(P_subdata)  # transform column 2-5 into a matrix
#rownames(mat_data) <- rnames 
head(mat_data)
single_column <- gather(subdata,key = GO, value = foldchanges,na.rm =T)
head(single_column)
quantile(single_column$foldchanges,probs = seq(0,1,0.01))
#0%        1%        2%        3%        4%        5%        6%        7%        8%        9%       10% 
#1.303440  1.304030  1.304620  1.305210  1.320673  1.338261  1.355848  1.360245  1.360245  1.360245  1.372196 
#11%       12%       13%       14%       15%       16%       17%       18%       19%       20%       21% 
#1.391318  1.410439  1.420000  1.420000  1.420000  1.448800  1.525600  1.602400  1.660000  1.660000  1.660000 
#22%       23%       24%       25%       26%       27%       28%       29%       30%       31%       32% 
#1.668800  1.739200  1.809600  1.880000  1.880000  1.880000  1.880000  1.891200  1.904000  1.916800  1.970400 
#33%       34%       35%       36%       37%       38%       39%       40%       41%       42%       43% 
#2.037600  2.104800  2.452000  2.967200  3.482400  3.794400  3.903200  4.012000  4.088007  4.109360  4.130713 
#44%       45%       46%       47%       48%       49%       50%       51%       52%       53%       54% 
#4.173389  4.280036  4.386684  4.480800  4.487200  4.493600  4.500000  4.615200  4.730400  4.845600  4.882400 
#55%       56%       57%       58%       59%       60%       61%       62%       63%       64%       65% 
#4.908000  4.933600  5.026400  5.141600  5.256800  5.516000  5.861600  6.207200  6.439200  6.557600  6.676000 
#66%       67%       68%       69%       70%       71%       72%       73%       74%       75%       76% 
#6.753600  6.763200  6.772800  6.807251  6.916257  7.025263  7.141817  7.311212  7.480606  7.650000  7.650000 
#77%       78%       79%       80%       81%       82%       83%       84%       85%       86%       87% 
#7.650000  7.650000  7.669600  7.692000  7.714400  7.787200  7.876800  7.966400  8.236000  8.613600  8.991200 
#88%       89%       90%       91%       92%       93%       94%       95%       96%       97%       98% 
#9.311200  9.573600  9.836000 10.025200 10.092400 10.159600 10.413200 11.226000 12.038800 12.969200 14.722800 
#99%      100% 
#  16.476400 18.230000 

#down
#0%      1%      2%      3%      4%      5%      6%      7%      8%      9%     10%     11%     12%     13% 
#  1.1000  1.1500  1.1600  1.1651  1.1768  1.1800  1.1900  1.2119  1.2336  1.2600  1.2600  1.2700  1.2800  1.2900 
#14%     15%     16%     17%     18%     19%     20%     21%     22%     23%     24%     25%     26%     27% 
#  1.2900  1.3000  1.3000  1.3000  1.3106  1.3200  1.3240  1.3500  1.3500  1.3782  1.3908  1.4000  1.4100  1.4159 
#28%     29%     30%     31%     32%     33%     34%     35%     36%     37%     38%     39%     40%     41% 
#  1.4200  1.4293  1.4400  1.4500  1.4644  1.4861  1.4978  1.5100  1.5200  1.5358  1.5700  1.5863  1.6100  1.6100 
#42%     43%     44%     45%     46%     47%     48%     49%     50%     51%     52%     53%     54%     55% 
#  1.6314  1.6431  1.6548  1.6700  1.6882  1.7099  1.7300  1.7400  1.7500  1.7667  1.7884  1.8802  1.9126  1.9870 
#56%     57%     58%     59%     60%     61%     62%     63%     64%     65%     66%     67%     68%     69% 
#  2.0152  2.0369  2.1786  2.1900  2.2300  2.2948  2.3254  2.3500  2.3964  2.4305  2.4754  2.5300  2.5412  2.5846 
#70%     71%     72%     73%     74%     75%     76%     77%     78%     79%     80%     81%     82%     83% 
#  2.6200  2.6514  2.7324  2.8441  2.8616  2.9450  3.0952  3.1627  3.2000  3.2644  3.3800  3.4000  3.5552  3.7411 
#84%     85%     86%     87%     88%     89%     90%     91%     92%     93%     94%     95%     96%     97% 
#  3.7824  3.8400  4.0386  4.1411  4.4204  4.6069  4.9430  5.0300  5.0300  5.1200  5.1700  5.6170  6.8796  7.2400 
#98%     99%    100% 
#7.5962 10.2327 12.0700 

#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.48,length=100),  # for red
#               seq(1.49,2.47,length=100),           # for yellow
#               seq(2.48,12.07,length=100))             # for green
#my_palette <- colorRampPalette(c("white","yellow", "red"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.47,length=100),  # for red
#               seq(1.48,1.90,length=100),           # for yellow
#               seq(1.91,10.7,length=100))             # for green
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(0,1.47,length=100),  # for red
#               seq(1.48,1.90,length=100),           # for yellow
#               seq(1.91,10.7,length=100))             # for green
#col_breaks = c(seq(0,2.34,length=100),  # for red
#               seq(2.35,3.73,length=100),           # for yellow
#               seq(3.74,14.22,length=100)) 
my_palette <- colorRampPalette(c("white","yellow", "red"))(n = 299)
#my_palette <- colorRampPalette(c("white","green", "blue"))(n = 299)
col_breaks = c(seq(0,1.45,length=100),  # for red
               seq(1.46,1.97,length=100),           # for yellow
               seq(1.98,7.52,length=100)) 
# creates a 5 x 5 inch image
pdf("~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/AgriGO/heatmaps_drought_up_P.pdf",    # create PNG for the heat map        
    width = 10,        # 5 x 300 pixels
    height =15
)        # smaller font size
dist_no_na <- function(mat) {
  edist <- dist(mat)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}
heatmap.2(mat_data,
          cellnote = mat_data,  # same data set for cell labels
          main = "Foldchange", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,6),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="none",# Don't draw dendrogram
          distfun=dist_no_na,
          cexRow=2,
          cexCol=2,
          Colv="NA",
          key =T,
          lwid = c(0.6, 4, 0.6),
          lmat = rbind( c(0, 3, 0), c(2, 1, 0), c(0, 4, 0) ),
          lhei = c(0.43, 2.6, 0.6))            # turn off column clustering

dev.off()               # close the PNG device


