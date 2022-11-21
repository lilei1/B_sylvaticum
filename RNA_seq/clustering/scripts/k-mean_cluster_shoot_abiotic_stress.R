library(edgeR)
library(SummarizedExperiment)


data.abiotic <- read.delim(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/normalized_VST_final_shoot_abiotic_formal.txt", sep="\t",row.names=1, header=T)
#data.abiotic <- read.delim(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/combined/normalized_VST_final_shoot_abiotic_freeze_tissue_novargen.txt", sep="\t",row.names=1, header=T)
head(data.abiotic)
#data.abiotic <- data.abiotic[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
#                                30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,
#                                16,17,18,19,20,21,22,23,24,25,26,27,28,29)]
#rownames(data.abiotic) <- data.abiotic$GeneID
#shoot <- data.abiotic[,c(1:59)]
shoot <- data.abiotic[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                                30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,
                                16,17,18,19,20,21,22,23,24,25,26,27,28,29)]

DEGs_shoot <- read.delim(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/all_sig_DEGs_across_stresses_time_gene.list.txt", sep="\t",row.names=1, header=T)
head(DEGs_shoot)
row.num <- which( rownames(shoot) %in% DEGs_shoot$x)
shoot <- shoot[row.num,]
nrow(shoot)
#[1] 17184
head(shoot)

shoot_var <- apply(shoot, 1, var)
quantile(shoot_var,probs=seq(0,1,0.01))
#0%          1%          2%          3%          4%          5%          6%          7% 
#  0.07241415  0.11560577  0.12787253  0.13961117  0.14865254  0.15655186  0.16407195  0.17186396 
#8%          9%         10%         11%         12%         13%         14%         15% 
#  0.17811227  0.18484667  0.19235254  0.19872775  0.20644275  0.21274184  0.21972487  0.22594674 
#16%         17%         18%         19%         20%         21%         22%         23% 
#  0.23160175  0.23881335  0.24510211  0.25211564  0.25874910  0.26576537  0.27237654  0.27924182 
#24%         25%         26%         27%         28%         29%         30%         31% 
#  0.28551911  0.29194723  0.29879480  0.30604904  0.31300846  0.32256770  0.32993495  0.33689110 
#32%         33%         34%         35%         36%         37%         38%         39% 
#  0.34419537  0.35171785  0.35908651  0.36610768  0.37430651  0.38187833  0.39051392  0.39939257 
#40%         41%         42%         43%         44%         45%         46%         47% 
#  0.40780654  0.41628473  0.42564352  0.43422895  0.44299709  0.45329149  0.46244474  0.47169534 
#48%         49%         50%         51%         52%         53%         54%         55% 
#  0.48186037  0.49139360  0.50220329  0.51316527  0.52496780  0.53780998  0.54930324  0.56115837 
#56%         57%         58%         59%         60%         61%         62%         63% 
#  0.57354682  0.58629876  0.60197701  0.61642131  0.63166156  0.65050173  0.66661507  0.68257676 
#64%         65%         66%         67%         68%         69%         70%         71% 
#  0.69978707  0.71656203  0.73844841  0.75823297  0.77909787  0.80172319  0.82566152  0.85090453 
#72%         73%         74%         75%         76%         77%         78%         79% 
#  0.87493825  0.90202762  0.92879390  0.95598009  0.98994818  1.02151300  1.05975331  1.09110776 
#80%         81%         82%         83%         84%         85%         86%         87% 
#  1.13548698  1.17563956  1.22339036  1.27451661  1.32739970  1.39631844  1.46948706  1.54202953 
#88%         89%         90%         91%         92%         93%         94%         95% 
#  1.62700596  1.72928661  1.83179230  1.96901066  2.10927280  2.30602338  2.51827932  2.78843911 
#96%         97%         98%         99%        100% 
#3.22701928  3.79699322  4.59798648  6.40121628 28.77072456 

shoot_mean <- apply(shoot, 1, mean)
quantile(shoot_mean,probs=seq(0,1,0.01))
#0%        1%        2%        3%        4%        5%        6%        7%        8%        9% 
#  2.430442  2.817136  3.055151  3.228100  3.376130  3.474811  3.559908  3.645344  3.720398  3.789247 
#10%       11%       12%       13%       14%       15%       16%       17%       18%       19% 
#  3.862313  3.924755  3.984017  4.040221  4.099003  4.162663  4.231690  4.303630  4.377962  4.447441 
#20%       21%       22%       23%       24%       25%       26%       27%       28%       29% 
#  4.512814  4.580221  4.650223  4.713488  4.789985  4.858663  4.935794  5.020184  5.083177  5.157268 
#30%       31%       32%       33%       34%       35%       36%       37%       38%       39% 
#  5.230944  5.296705  5.363778  5.431108  5.492579  5.567126  5.642487  5.721243  5.792276  5.856979 
#40%       41%       42%       43%       44%       45%       46%       47%       48%       49% 
#  5.926794  6.001715  6.072432  6.139156  6.211481  6.286411  6.361009  6.441061  6.525234  6.591706 
#50%       51%       52%       53%       54%       55%       56%       57%       58%       59% 
#  6.670000  6.737738  6.818435  6.893601  6.969154  7.031550  7.094955  7.168641  7.242811  7.322633 
#60%       61%       62%       63%       64%       65%       66%       67%       68%       69% 
#  7.397893  7.468854  7.543184  7.625962  7.702741  7.782899  7.855998  7.931547  8.011384  8.088687 
#70%       71%       72%       73%       74%       75%       76%       77%       78%       79% 
#  8.166962  8.257742  8.335193  8.408803  8.485796  8.570780  8.647548  8.731458  8.812110  8.905613 
#80%       81%       82%       83%       84%       85%       86%       87%       88%       89% 
#  8.991739  9.081065  9.181152  9.283730  9.382790  9.476053  9.583129  9.696136  9.807345  9.931286 
#90%       91%       92%       93%       94%       95%       96%       97%       98%       99% 
#  10.068270 10.193099 10.354287 10.517206 10.735639 10.975547 11.249883 11.611038 12.106833 12.936832 
#100% 
#18.275886 
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/filter_abiotic.pdf", width = 15,height = 13)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/filter_shoot_abiotic.pdf", width = 15,height = 13)

plot(log2(shoot_mean), log2(shoot_var), pch='.',cex.axis=2)
abline(h=log2(1.5), col='red')#tap 26% 
abline(v=log2(4), col='red')#top 73%
text(x=3,y=4, font=30,labels="variance > 1.5 &\n mean > 4", col='red')
dev.off()

 
shoot <- shoot[which(shoot_var > 1.5 & shoot_mean > 4),]
nrow(shoot)
#[1] 1911
head(shoot) 

scaledata <- t(scale(t(shoot)))
head(scaledata)
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
TreeC = as.dendrogram(hc, method="average")
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/cluster_abiotic.pdf", width = 13,height = 17)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/cluster_shoot_abiotic.pdf", width = 13,height = 17)
plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")
dev.off()
#loadedNamespaces()

pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/SSE_shoot_abiotic.pdf", width = 10,height = 10)
#library(cluster)
wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
#wss <- vector('numeric', 20)

for (i in 2:20) wss[i] <- sum(kmeans(scaledata,centers=i)$withinss)

plot(1:20, wss, type="b",cex.axis=2,lwd=2,xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(v=5, col='red')#top 73%
dev.off()

 
##cluster =7

library(cluster)
sil <- rep(0, 20)
#repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
  ss <- silhouette(k1to20$cluster, dist(scaledata))
  sil[i] <- mean(ss[, 3])
}
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/Av_silhouette_abiotic.pdf", width = 10,height = 10)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/Av_silhouette_shoot_abiotic.pdf", width = 10,height = 10)

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)
cat("Average silhouette width optimal number of clusters:", which.max(sil), "\n")
dev.off()
#2-3 cluster could be the same

BiocManager::install("vegan")
library(vegan)
fit <- cascadeKM(scaledata, 1, 20, iter = 100)
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/Calinsky_abiotic.pdf", width = 10,height = 10)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/Calinsky_shoot_abiotic.pdf", width = 10,height = 10)

plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")
dev.off()
##cluster =2

library(cluster)
set.seed(13)
gap <- clusGap(scaledata, kmeans, 20, B = 100, verbose = interactive())
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/Gap_abiotic.pdf", width = 10,height = 10)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/Gap_shoot_abiotic.pdf", width = 10,height = 10)

plot(gap, main = "Gap statistic")
abline(v=which.max(gap$Tab[,3]), lty = 2)
dev.off()
##cluster=20

BiocManager::install("apcluster")
library(apcluster)
d.apclus <- apcluster(negDistMat(r=2), scaledata)
cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")
## affinity propogation optimal number of clusters: 63
#uncomment the next line for the heatmap:
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/Affinity_abiotic.pdf", width = 15,height = 15)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/Affinity_shoot_abiotic.pdf", width = 15,height = 15)

heatmap(d.apclus,cexRow=0, cexCol=0)
dev.off()
#cluster=8 is good

BiocManager::install("gplots")
unloadNamespace("GeneOverlap")
library(gplots)

#make the matrix
dist <- cor(t(scaledata), method="pearson")
#make the tree
hr <- hclust(as.dist(1-dist), method="complete") # Cluster rows by Pearson correlation.
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/Hier_cluster_abiotic.pdf", width = 15,height = 15)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/Hier_shoot_abiotic.pdf", width = 15,height = 15)

heatmap.2(dist,
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hr),
          scale="row",
          margins = c(2, 2),
          cexCol = 0.7,
          #labRow = F,
          #labCol = F,
          main = "Heatmap",
          trace = "none"
)
dev.off()
#cluster =8

set.seed(20)
kClust <- kmeans(scaledata, centers=8, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)
library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/final_cluster5_abiotic.pdf", width = 20,height = 13)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/final_cluster8_shoot_abiotic.pdf", width = 20,height = 13)

#plot
#p1 <- ggplot(Kmolten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
#                                                    "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
#                                                    "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
#                                                    "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
#                                                    "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
#                                                    "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
#                                                    "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
#                                                    "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
#                                                    "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
#                                                    "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
#                                                    "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
#                                                    "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
#                                                    "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
#                                                    "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
#                                                    "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
#                                                    "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
#                                                    "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
#                                                    "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
#                                                    "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
#                                                    "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
#                                                  )),y=value, group=cluster, colour=as.factor(cluster))) + 
#  geom_point() + 
#  geom_line() +
#  xlab("Time") +
#  ylab("Expression") +
#  labs(title= "Cluster Expression by Time",color = "Cluster")
p1 <- ggplot(Kmolten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                    "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                    "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                    "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                    "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                    "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                    "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                    "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                    "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                    "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                    "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                    "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                    "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                    "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                    "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                    "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                    "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                    "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                    "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                    "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
                                                    )),y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Treatments") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")

p2 <- p1+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p2)
dev.off()
cor(kClustcentroids)
#K=8 for freeze tissue samples
#1          2          3           4           5           6           7           8
#1  1.0000000  0.4712493 -0.8593979  0.91166253 -0.74771100 -0.37062965 -0.87769294 -0.11830284
#2  0.4712493  1.0000000 -0.1060483  0.73356570 -0.91835755  0.63999595 -0.77811481  0.26164585
#3 -0.8593979 -0.1060483  1.0000000 -0.74988430  0.48680967  0.65973035  0.70582383  0.45133740
#4  0.9116625  0.7335657 -0.7498843  1.00000000 -0.93783168 -0.03379646 -0.99587210 -0.09773738
#5 -0.7477110 -0.9183575  0.4868097 -0.93783168  1.00000000 -0.30153321  0.96156119 -0.06163483
#6 -0.3706296  0.6399959  0.6597303 -0.03379646 -0.30153321  1.00000000 -0.03847582  0.40387819
#7 -0.8776929 -0.7781148  0.7058238 -0.99587210  0.96156119 -0.03847582  1.00000000  0.07958515
#8 -0.1183028  0.2616458  0.4513374 -0.09773738 -0.06163483  0.40387819  0.07958515  1.00000000

#Seems OK!!


#try cluster3
#Subset the cores molten dataframe so we can plot the core
core3 <- Kmolten[Kmolten$cluster=="3",]

#get cluster 1
K3 <- (scaledata[kClusters==3,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core3$value)}
score <- apply(K3, 1, corscore)
#get the data frame into long format for plotting
K3molten <- melt(K3)
colnames(K3molten) <- c('gene','sample','value')
#add the score
K3molten <- merge(K3molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K3molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K3molten$order_factor <- 1:length(K3molten$gene)
#order the dataframe by score
K3molten <- K3molten[order(K3molten$score),]
#set the order by setting the factors
K3molten$order_factor <- factor(K3molten$order_factor , levels = K3molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the3cluster_shoot_abiotic.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K3molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                       "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                       "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                       "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                       "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                       "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                       "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                       "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                       "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                       "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                       "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                       "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                       "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                       "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                       "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                       "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                       "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                       "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                       "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                       "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
#this adds the core 
  geom_line(data=core3, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 3 Expression by Time",color = "Score")
#p3 <- ggplot(K3molten, aes(x= factor (sample,level=c("leaf_cold_1","leaf_cold_2","leaf_cold_3",
#                                                     "leaf_freeze_1","leaf_freeze_2","leaf_freeze_3",
#                                                     "leaf_recovery_1","leaf_recovery_2","leaf_recovery_3",
#                                                     "crown_cold_1","crown_cold_2","crown_cold_3",
#                                                     "crown_freeze_1","crown_freeze_2","crown_freeze_3",
#                                                     "crown_recovery_1","crown_recovery_2","crown_recovery_3"
#)),y=value))+ 
#  geom_line(aes(colour=score, group=gene)) +
#  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
#  geom_line(data=core3, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
#  xlab("Time") +
#  ylab("Expression") +
#  labs(title= "Cluster 3 Expression by Time",color = "Score")



p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

##try cluster5
#Subset the cores molten dataframe so we can plot the core
core5 <- Kmolten[Kmolten$cluster=="5",]

#get cluster 5
K5 <- (scaledata[kClusters==5,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core5$value)}
score <- apply(K5, 1, corscore)
#get the data frame into long format for plotting
K5molten <- melt(K5)
colnames(K5molten) <- c('gene','sample','value')
#add the score
K5molten <- merge(K5molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K5molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K5molten$order_factor <- 1:length(K5molten$gene)
#order the dataframe by score
K5molten <- K5molten[order(K5molten$score),]
#set the order by setting the factors
K5molten$order_factor <- factor(K5molten$order_factor , levels = K5molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the5cluster_shoot_abiotic.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K5molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core5, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 5 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster1
#Subset the cores molten dataframe so we can plot the core
core1 <- Kmolten[Kmolten$cluster=="1",]

#get cluster 1
K1 <- (scaledata[kClusters==1,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core1$value)}
score <- apply(K1, 1, corscore)

#get the data frame into long format for plotting
K1molten <- melt(K1)
colnames(K1molten) <- c('gene','sample','value')
#add the score
K1molten <- merge(K1molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K1molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K1molten$order_factor <- 1:length(K1molten$gene)
#order the dataframe by score
K1molten <- K1molten[order(K1molten$score),]
#set the order by setting the factors
K1molten$order_factor <- factor(K1molten$order_factor , levels = K1molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the1cluster_shoot_abiotic.pdf", width = 20,height = 13)
# Everything on the same plot
p3 <- ggplot(K1molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core1, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 1 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()
#try cluster2
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster=="2",]

#get cluster 1
K2 <- (scaledata[kClusters==2,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the2cluster_shoot_abiotic.pdf", width = 20,height = 13)
p3 <- ggplot(K2molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core2, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 2 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster4
#Subset the cores molten dataframe so we can plot the core
core4 <- Kmolten[Kmolten$cluster=="4",]

#get cluster 4
K4 <- (scaledata[kClusters==4,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core4$value)}
score <- apply(K4, 1, corscore)
#get the data frame into long format for plotting
K4molten <- melt(K4)
colnames(K4molten) <- c('gene','sample','value')
#add the score
K4molten <- merge(K4molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K4molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K4molten$order_factor <- 1:length(K4molten$gene)
#order the dataframe by score
K4molten <- K4molten[order(K4molten$score),]
#set the order by setting the factors
K4molten$order_factor <- factor(K4molten$order_factor , levels = K4molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the4cluster_shoot_abiotic.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K4molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core4, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 4 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster6
#Subset the cores molten dataframe so we can plot the core
core6 <- Kmolten[Kmolten$cluster=="6",]
#get cluster 6
K6 <- (scaledata[kClusters==6,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core6$value)}
score <- apply(K6, 1, corscore)
#get the data frame into long format for plotting
K6molten <- melt(K6)
colnames(K6molten) <- c('gene','sample','value')
#add the score
K6molten <- merge(K6molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K6molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K6molten$order_factor <- 1:length(K6molten$gene)
#order the dataframe by score
K6molten <- K6molten[order(K6molten$score),]
#set the order by setting the factors
K6molten$order_factor <- factor(K6molten$order_factor , levels = K6molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the6cluster_shoot_abiotic.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K6molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core6, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 6 Expression by Time",color = "Score")

p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster7
#Subset the cores molten dataframe so we can plot the core
core7 <- Kmolten[Kmolten$cluster=="7",]
#get cluster 7
K7 <- (scaledata[kClusters==7,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core7$value)}
score <- apply(K7, 1, corscore)
#get the data frame into long format for plotting
K7molten <- melt(K7)
colnames(K7molten) <- c('gene','sample','value')
#add the score
K7molten <- merge(K7molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K7molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K7molten$order_factor <- 1:length(K7molten$gene)
#order the dataframe by score
K7molten <- K7molten[order(K7molten$score),]
#set the order by setting the factors
K7molten$order_factor <- factor(K7molten$order_factor , levels = K7molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the7cluster_shoot_abiotic.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K7molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core7, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 7 Expression by Time",color = "Score")

p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster8
#Subset the cores molten dataframe so we can plot the core
core8 <- Kmolten[Kmolten$cluster=="8",]
#get cluster 6
K8 <- (scaledata[kClusters==8,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core8$value)}
score <- apply(K8, 1, corscore)
#get the data frame into long format for plotting
K8molten <- melt(K8)
colnames(K8molten) <- c('gene','sample','value')
#add the score
K8molten <- merge(K8molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K8molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K8molten$order_factor <- 1:length(K8molten$gene)
#order the dataframe by score
K8molten <- K8molten[order(K8molten$score),]
#set the order by setting the factors
K8molten$order_factor <- factor(K8molten$order_factor , levels = K8molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the8cluster_shoot_abiotic.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K8molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core8, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 8 Expression by Time",color = "Score")

p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

##Since cluster 1,3,5 did not good, so I have to redo it with 5 clusters!!!

set.seed(20)
kClust <- kmeans(scaledata, centers=5, nstart = 1000, iter.max = 20)
kClusters <- kClust$cluster
# function to find centroid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters)
library(ggplot2)
library(reshape)
#get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c('sample','cluster','value')
#pdf(file="~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/timecourse/final_cluster5_abiotic.pdf", width = 20,height = 13)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/final_cluster5_shoot_abiotic.pdf", width = 20,height = 13)

#plot
#p1 <- ggplot(Kmolten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
#                                                    "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
#                                                    "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
#                                                    "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
#                                                    "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
#                                                    "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
#                                                    "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
#                                                    "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
#                                                    "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
#                                                    "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
#                                                    "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
#                                                    "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
#                                                    "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
#                                                    "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
#                                                    "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
#                                                    "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
#                                                    "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
#                                                    "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
#                                                    "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
#                                                    "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
#                                                  )),y=value, group=cluster, colour=as.factor(cluster))) + 
#  geom_point() + 
#  geom_line() +
#  xlab("Time") +
#  ylab("Expression") +
#  labs(title= "Cluster Expression by Time",color = "Cluster")
p1 <- ggplot(Kmolten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                    "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                    "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                    "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                    "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                    "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                    "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                    "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                    "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                    "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                    "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                    "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                    "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                    "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                    "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                    "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                    "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                    "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                    "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                    "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value, group=cluster, colour=as.factor(cluster))) + 
  geom_point() + 
  geom_line() +
  xlab("Treatments") +
  ylab("Expression") +
  labs(title= "Cluster Expression by Time",color = "Cluster")

p2 <- p1+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p2)
dev.off()
cor(kClustcentroids)
#1          2          3          4          5
#1  1.0000000  0.1378149 -0.2615647 -0.7200599 -0.1896760
#2  0.1378149  1.0000000 -0.9645862 -0.1191489  0.8111491
#3 -0.2615647 -0.9645862  1.0000000  0.1829027 -0.7945225
#4 -0.7200599 -0.1191489  0.1829027  1.0000000  0.2838667
#5 -0.1896760  0.8111491 -0.7945225  0.2838667  1.0000000
#Seems OK!!


#try cluster3
#Subset the cores molten dataframe so we can plot the core
core3 <- Kmolten[Kmolten$cluster=="3",]

#get cluster 1
K3 <- (scaledata[kClusters==3,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core3$value)}
score <- apply(K3, 1, corscore)
#get the data frame into long format for plotting
K3molten <- melt(K3)
colnames(K3molten) <- c('gene','sample','value')
#add the score
K3molten <- merge(K3molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K3molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K3molten$order_factor <- 1:length(K3molten$gene)
#order the dataframe by score
K3molten <- K3molten[order(K3molten$score),]
#set the order by setting the factors
K3molten$order_factor <- factor(K3molten$order_factor , levels = K3molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the3cluster_shoot_abiotic_all5cluster.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K3molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core3, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 3 Expression by Time",color = "Score")
#p3 <- ggplot(K3molten, aes(x= factor (sample,level=c("leaf_cold_1","leaf_cold_2","leaf_cold_3",
#                                                     "leaf_freeze_1","leaf_freeze_2","leaf_freeze_3",
#                                                     "leaf_recovery_1","leaf_recovery_2","leaf_recovery_3",
#                                                     "crown_cold_1","crown_cold_2","crown_cold_3",
#                                                     "crown_freeze_1","crown_freeze_2","crown_freeze_3",
#                                                     "crown_recovery_1","crown_recovery_2","crown_recovery_3"
#)),y=value))+ 
#  geom_line(aes(colour=score, group=gene)) +
#  scale_colour_gradientn(colours=c('blue1','red2')) +
#this adds the core 
#  geom_line(data=core3, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
#  xlab("Time") +
#  ylab("Expression") +
#  labs(title= "Cluster 3 Expression by Time",color = "Score")



p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

##try cluster5
#Subset the cores molten dataframe so we can plot the core
core5 <- Kmolten[Kmolten$cluster=="5",]

#get cluster 5
K5 <- (scaledata[kClusters==5,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core5$value)}
score <- apply(K5, 1, corscore)
#get the data frame into long format for plotting
K5molten <- melt(K5)
colnames(K5molten) <- c('gene','sample','value')
#add the score
K5molten <- merge(K5molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K5molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K5molten$order_factor <- 1:length(K5molten$gene)
#order the dataframe by score
K5molten <- K5molten[order(K5molten$score),]
#set the order by setting the factors
K5molten$order_factor <- factor(K5molten$order_factor , levels = K5molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the5cluster_shoot_abiotic_all5cluster.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K5molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core5, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 5 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster1
#Subset the cores molten dataframe so we can plot the core
core1 <- Kmolten[Kmolten$cluster=="1",]

#get cluster 1
K1 <- (scaledata[kClusters==1,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core1$value)}
score <- apply(K1, 1, corscore)

#get the data frame into long format for plotting
K1molten <- melt(K1)
colnames(K1molten) <- c('gene','sample','value')
#add the score
K1molten <- merge(K1molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K1molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K1molten$order_factor <- 1:length(K1molten$gene)
#order the dataframe by score
K1molten <- K1molten[order(K1molten$score),]
#set the order by setting the factors
K1molten$order_factor <- factor(K1molten$order_factor , levels = K1molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the1cluster_shoot_abiotic_all5cluster.pdf", width = 20,height = 13)
# Everything on the same plot
p3 <- ggplot(K1molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core1, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 1 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()
#try cluster2
#Subset the cores molten dataframe so we can plot the core
core2 <- Kmolten[Kmolten$cluster=="2",]

#get cluster 1
K2 <- (scaledata[kClusters==2,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)
#get the data frame into long format for plotting
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
#add the score
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K2molten$order_factor <- 1:length(K2molten$gene)
#order the dataframe by score
K2molten <- K2molten[order(K2molten$score),]
#set the order by setting the factors
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the2cluster_shoot_abiotic_all5cluster.pdf", width = 20,height = 13)
p3 <- ggplot(K2molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core2, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 2 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

#try cluster4
#Subset the cores molten dataframe so we can plot the core
core4 <- Kmolten[Kmolten$cluster=="4",]

#get cluster 4
K4 <- (scaledata[kClusters==4,])
#calculate the correlation with the core
corscore <- function(x){cor(x,core4$value)}
score <- apply(K4, 1, corscore)
#get the data frame into long format for plotting
K4molten <- melt(K4)
colnames(K4molten) <- c('gene','sample','value')
#add the score
K4molten <- merge(K4molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K4molten) <- c('gene','sample','value','score')
#order the dataframe by score
#to do this first create an ordering factor
K4molten$order_factor <- 1:length(K4molten$gene)
#order the dataframe by score
K4molten <- K4molten[order(K4molten$score),]
#set the order by setting the factors
K4molten$order_factor <- factor(K4molten$order_factor , levels = K4molten$order_factor)
pdf(file="~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/the4cluster_shoot_abiotic_all5cluster.pdf", width = 20,height = 13)

# Everything on the same plot
p3 <- ggplot(K4molten, aes(x= factor (sample,level=c("Shoot_CK_1h.s1","Shoot_CK_1h.s2","Shoot_CK_1h.s3",
                                                     "Shoot_CK_2h.s1","Shoot_CK_2h.s2","Shoot_CK_2h.s3",
                                                     "Shoot_CK_5h.s1","Shoot_CK_5h.s2","Shoot_CK_5h.s3",
                                                     "Shoot_CK_10h.s1","Shoot_CK_10h.s2","Shoot_CK_10h.s3",
                                                     "Shoot_CK_24h.s1","Shoot_CK_24h.s2","Shoot_CK_24h.s3",
                                                     "Shoot_drought_1h.s1","Shoot_drought_1h.s2","Shoot_drought_1h.s3",
                                                     "Shoot_drought_2h.s1","Shoot_drought_2h.s2","Shoot_drought_2h.s3",
                                                     "Shoot_drought_5h.s1","Shoot_drought_5h.s2","Shoot_drought_5h.s3",
                                                     "Shoot_drought_10h.s1","Shoot_drought_10h.s2","Shoot_drought_10h.s3",
                                                     "Shoot_drought_24h.s1","Shoot_drought_24h.s2","Shoot_drought_24h.s3",
                                                     "Shoot_heat_1h.s1","Shoot_heat_1h.s2","Shoot_heat_1h.s3",
                                                     "Shoot_heat_2h.s1","Shoot_heat_2h.s2","Shoot_heat_2h.s3",
                                                     "Shoot_heat_5h.s1","Shoot_heat_5h.s2","Shoot_heat_5h.s3",
                                                     "Shoot_heat_10h.s1","Shoot_heat_10h.s2","Shoot_heat_10h.s3",
                                                     "Shoot_heat_24h.s1","Shoot_heat_24h.s2","Shoot_heat_24h.s3",
                                                     "Shoot_salt_1h.s1","Shoot_salt_1h.s2","Shoot_salt_1h.s3",
                                                     "Shoot_salt_2h.s1","Shoot_salt_2h.s2","Shoot_salt_2h.s3",
                                                     "Shoot_salt_5h.s1","Shoot_salt_5h.s2","Shoot_salt_5h.s3",
                                                     "Shoot_salt_10h.s1","Shoot_salt_10h.s2","Shoot_salt_10h.s3",
                                                     "Shoot_salt_24h.s1","Shoot_salt_24h.s2"
)),y=value)) + 
  geom_line(aes(colour=score, group=gene)) +
  scale_colour_gradientn(colours=c('blue1','red2')) +
  #this adds the core 
  geom_line(data=core4, aes(sample,value, group=cluster), color="black",inherit.aes=FALSE) +
  xlab("Time") +
  ylab("Expression") +
  labs(title= "Cluster 4 Expression by Time",color = "Score")


p4 <- p3+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p4)
dev.off()

head(K3molten)
write.table(K3molten, file = "~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/K3molten_cluster_shoot_abiotic.txt", 
            sep="\t", quote = FALSE)
write.table(K4molten, file = "~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/K4molten_cluster_shoot_abiotic.txt", 
            sep="\t", quote = FALSE)
write.table(K5molten, file = "~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/K5molten_cluster_shoot_abiotic.txt", 
            sep="\t", quote = FALSE)
write.table(K1molten, file = "~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/K1molten_cluster_shoot_abiotic.txt", 
            sep="\t", quote = FALSE)
write.table(K2molten, file = "~/MorrellLab Dropbox/Li Lei/Brachpodium/sylvaticum/RNA_seq/timecourse/K2molten_cluster_shoot_abiotic.txt", 
            sep="\t", quote = FALSE)


#write.table(K6molten, file = "~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/combined/K6molten_cluster_freeze.txt", 
            sep="\t", quote = FALSE)
#write.table(K7molten, file = "~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/combined/K7molten_cluster_freeze.txt", 
            sep="\t", quote = FALSE)
#write.table(K8molten, file = "~/Dropbox (Personal)/Brachpodium/sylvaticum/RNA_seq/combined/K8molten_cluster_freeze.txt", 
            sep="\t", quote = FALSE)
