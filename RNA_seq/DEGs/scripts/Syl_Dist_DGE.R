#Read the comparable data from Distachyon and sylvaticum:

all <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/Syl-Dist_matched_DGE.txt",header=T)
head(all)
all$distID
nrow(all)
#heat1h
heat1h <- all[,c("Syl_ID","S_heat1h.log2FoldChange","S_heat1h.padj","distID","D_heat1h.log2FoldChange","D_heat1h.padj")]
head(heat1h)
nrow(heat1h)
#[1] 20735
heat1h <- na.omit(heat1h)
nrow(heat1h)
#16231
heat1h_both_DGE <- heat1h[(heat1h$S_heat1h.padj<0.05 & heat1h$D_heat1h.padj<0.05),]
nrow(heat1h_both_DGE)
#3233
#heat1h_both_DGE$distID
heat1h_both_DGE_high <- heat1h_both_DGE[(abs(heat1h_both_DGE$S_heat1h.log2FoldChange)>=4 & abs(heat1h_both_DGE$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_both_DGE_high)
#16
#both negative
heat1h_both_negative <- heat1h_both_DGE[(heat1h_both_DGE$S_heat1h.log2FoldChange<0 & heat1h_both_DGE$D_heat1h.log2FoldChange<0),]
nrow(heat1h_both_negative)
#782
heat1h_both_negative_high <- heat1h_both_negative[(abs(heat1h_both_negative$S_heat1h.log2FoldChange)>=4 & abs(heat1h_both_negative$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_both_negative_high)
#0

#both positive
heat1h_both_positive <- heat1h_both_DGE[(heat1h_both_DGE$S_heat1h.log2FoldChange>0 & heat1h_both_DGE$D_heat1h.log2FoldChange>0),]
nrow(heat1h_both_positive)
#811
heat1h_both_positive_high <- heat1h_both_positive[(abs(heat1h_both_positive$S_heat1h.log2FoldChange)>=4 & abs(heat1h_both_positive$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_both_positive_high)
#8

#negative in S but positive in D
heat1h_pS_nD <- heat1h_both_DGE[(heat1h_both_DGE$S_heat1h.log2FoldChange>0 & heat1h_both_DGE$D_heat1h.log2FoldChange<0),]
nrow(heat1h_pS_nD)
#841
heat1h_pS_nD_high <- heat1h_pS_nD[(abs(heat1h_pS_nD$S_heat1h.log2FoldChange)>=4 & abs(heat1h_pS_nD$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_pS_nD_high)
#5
heat1h_nS_pD <- heat1h_both_DGE[(heat1h_both_DGE$S_heat1h.log2FoldChange<0 & heat1h_both_DGE$D_heat1h.log2FoldChange>0),]
nrow(heat1h_nS_pD)
#799

heat1h_nS_pD_high <- heat1h_nS_pD[(abs(heat1h_nS_pD$S_heat1h.log2FoldChange)>=4 & abs(heat1h_nS_pD$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_nS_pD_high)
#3

#look for the same direct and same fold changes
ratio <- heat1h_both_DGE$S_heat1h.log2FoldChange/heat1h_both_DGE$D_heat1h.log2FoldChange
heat1h_S_D_plasticity <- heat1h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(heat1h_S_D_plasticity)
#221 check the GO for those genes
heat1h_S_D_plasticity_high <- heat1h_S_D_plasticity[(abs(heat1h_S_D_plasticity$S_heat1h.log2FoldChange)>=4 & abs(heat1h_S_D_plasticity$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_S_D_plasticity_high)
#4

#make plot for the heat1h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat1h_DEG_S_D.pdf",width = 8,height = 8)
plot(heat1h_both_negative$S_heat1h.log2FoldChange,heat1h_both_negative$D_heat1h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-1h")
points(heat1h_both_positive$S_heat1h.log2FoldChange,heat1h_both_positive$D_heat1h.log2FoldChange,col="green")
points(heat1h_pS_nD$S_heat1h.log2FoldChange,heat1h_pS_nD$D_heat1h.log2FoldChange,col="blue")
points(heat1h_nS_pD$S_heat1h.log2FoldChange,heat1h_nS_pD$D_heat1h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat1h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(heat1h_both_negative_high$S_heat1h.log2FoldChange,heat1h_both_negative_high$D_heat1h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-1h")
points(heat1h_both_positive_high$S_heat1h.log2FoldChange,heat1h_both_positive_high$D_heat1h.log2FoldChange,col="green")
points(heat1h_pS_nD_high$S_heat1h.log2FoldChange,heat1h_pS_nD_high$D_heat1h.log2FoldChange,col="blue")
points(heat1h_nS_pD_high$S_heat1h.log2FoldChange,heat1h_nS_pD_high$D_heat1h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


heat1h_S_unique <- heat1h[(heat1h$S_heat1h.padj<0.05 & heat1h$D_heat1h.padj>0.05),]
nrow(heat1h_S_unique)
#4336
heat1h_S_unique_high <- heat1h_S_unique[(abs(heat1h_S_unique$S_heat1h.log2FoldChange)>=4),]
nrow(heat1h_S_unique_high)
#282
heat1h_D_unique <- heat1h[(heat1h$S_heat1h.padj>0.05 & heat1h$D_heat1h.padj<0.05),]
nrow(heat1h_D_unique)
#3837
heat1h_D_unique$distID
head(heat1h_D_unique)
heat1h_D_unique_high <- heat1h_D_unique[(abs(heat1h_D_unique$D_heat1h.log2FoldChange)>=4),]
nrow(heat1h_D_unique_high)
head(heat1h_D_unique_high)
#294
#head(heat1h_D_unique_high)
#heat1h_D_unique_high$distID
heat1h_both_NOT_DGE <- heat1h[(heat1h$S_heat1h.padj>0.05 & heat1h$D_heat1h.padj>0.05),]
nrow(heat1h_both_NOT_DGE)
#4825

#Heat2h
heat2h <- all[,c("Syl_ID","S_heat2h.log2FoldChange","S_heat2h.padj","distID","D_heat2h.log2FoldChange","D_heat2h.padj")]
head(heat2h)
nrow(heat2h)
#[1] 20735
heat2h <- na.omit(heat2h)
nrow(heat2h)
#15639
heat2h_both_DGE <- heat2h[(heat2h$S_heat2h.padj<0.05 & heat2h$D_heat2h.padj<0.05),]
nrow(heat2h_both_DGE)
#3397


#both negative
heat2h_both_negative <- heat2h_both_DGE[(heat2h_both_DGE$S_heat2h.log2FoldChange<0 & heat2h_both_DGE$D_heat2h.log2FoldChange<0),]
nrow(heat2h_both_negative)
#883
heat2h_both_negative_high <- heat2h_both_negative[(abs(heat2h_both_negative$S_heat2h.log2FoldChange)>=4 & abs(heat2h_both_negative$D_heat2h.log2FoldChange)>=4),]
nrow(heat2h_both_negative_high)
#0
#heat2h_both_negative_high <- heat2h_both_negative[(abs(heat2h_both_negative$D_heat2h.log2FoldChange)>=4),]

#both positive
heat2h_both_positive <- heat2h_both_DGE[(heat2h_both_DGE$S_heat2h.log2FoldChange>0 & heat2h_both_DGE$D_heat2h.log2FoldChange>0),]
nrow(heat2h_both_positive)
#813
heat2h_both_positive_high <- heat2h_both_positive[(abs(heat2h_both_positive$S_heat2h.log2FoldChange)>=4 & abs(heat2h_both_positive$D_heat2h.log2FoldChange)>=4),]
nrow(heat2h_both_positive_high)
#2

#negative in S but positive in D
heat2h_pS_nD <- heat2h_both_DGE[(heat2h_both_DGE$S_heat2h.log2FoldChange>0 & heat2h_both_DGE$D_heat2h.log2FoldChange<0),]
nrow(heat2h_pS_nD)
#839
heat2h_pS_nD_high <- heat2h_pS_nD[(abs(heat2h_pS_nD$S_heat2h.log2FoldChange)>=4 & abs(heat2h_pS_nD$D_heat2h.log2FoldChange)>=4),]
nrow(heat2h_pS_nD_high)
#3
heat2h_nS_pD <- heat2h_both_DGE[(heat2h_both_DGE$S_heat2h.log2FoldChange<0 & heat2h_both_DGE$D_heat2h.log2FoldChange>0),]
nrow(heat2h_nS_pD)
#862

heat2h_nS_pD_high <- heat2h_nS_pD[(abs(heat2h_nS_pD$S_heat2h.log2FoldChange)>=4 & abs(heat2h_nS_pD$D_heat2h.log2FoldChange)>=4),]
nrow(heat2h_nS_pD_high)
#3

#look for the same direct and same fold changes
ratio <- heat2h_both_DGE$S_heat2h.log2FoldChange/heat2h_both_DGE$D_heat2h.log2FoldChange
heat2h_S_D_plasticity <- heat2h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(heat2h_S_D_plasticity)
#265 check the GO for those genes
heat2h_S_D_plasticity_high <- heat2h_S_D_plasticity[(abs(heat2h_S_D_plasticity$S_heat2h.log2FoldChange)>=4 & abs(heat2h_S_D_plasticity$D_heat2h.log2FoldChange)>=4),]
nrow(heat2h_S_D_plasticity_high)
#1

#make plot for the heat2h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat2h_DEG_S_D.pdf",width = 8,height = 8)
plot(heat2h_both_negative$S_heat2h.log2FoldChange,heat2h_both_negative$D_heat2h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-2h")
points(heat2h_both_positive$S_heat2h.log2FoldChange,heat2h_both_positive$D_heat2h.log2FoldChange,col="green")
points(heat2h_pS_nD$S_heat2h.log2FoldChange,heat2h_pS_nD$D_heat2h.log2FoldChange,col="blue")
points(heat2h_nS_pD$S_heat2h.log2FoldChange,heat2h_nS_pD$D_heat2h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat2h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(heat2h_both_negative_high$S_heat2h.log2FoldChange,heat2h_both_negative_high$D_heat2h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-2h")
points(heat2h_both_positive_high$S_heat2h.log2FoldChange,heat2h_both_positive_high$D_heat2h.log2FoldChange,col="green")
points(heat2h_pS_nD_high$S_heat2h.log2FoldChange,heat2h_pS_nD_high$D_heat2h.log2FoldChange,col="blue")
points(heat2h_nS_pD_high$S_heat2h.log2FoldChange,heat2h_nS_pD_high$D_heat2h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#3233


heat2h_S_unique <- heat2h[(heat2h$S_heat2h.padj<0.05 & heat2h$D_heat2h.padj>0.05),]
nrow(heat2h_S_unique)
#4227
heat2h_S_unique_high <- heat2h_S_unique[(abs(heat2h_S_unique$S_heat2h.log2FoldChange)>=4),]
nrow(heat2h_S_unique_high)
head(heat2h_S_unique_high)
#180
heat2h_D_unique <- heat2h[(heat2h$S_heat2h.padj>0.05 & heat2h$D_heat2h.padj<0.05),]
nrow(heat2h_D_unique)
#3623
heat2h_D_unique_high <- heat2h_D_unique[(abs(heat2h_D_unique$D_heat2h.log2FoldChange)>=4),]
nrow(heat2h_D_unique_high)
#191
head(heat2h_D_unique_high)

heat2h_both_NOT_DGE <- heat2h[(heat2h$S_heat2h.padj>0.05 & heat2h$D_heat2h.padj>0.05),]
nrow(heat2h_both_NOT_DGE)
#4392

#5h:
heat5h <- all[,c("Syl_ID","S_heat5h.log2FoldChange","S_heat5h.padj","distID","D_heat5h.log2FoldChange","D_heat5h.padj")]
head(heat5h)
nrow(heat5h)
#[1] 20735
heat5h <- na.omit(heat5h)
nrow(heat5h)
#15013
heat5h_both_DGE <- heat5h[(heat5h$S_heat5h.padj<0.05 & heat5h$D_heat5h.padj<0.05),]
nrow(heat5h_both_DGE)
#2017


#both negative
heat5h_both_negative <- heat5h_both_DGE[(heat5h_both_DGE$S_heat5h.log2FoldChange<0 & heat5h_both_DGE$D_heat5h.log2FoldChange<0),]
nrow(heat5h_both_negative)
#449
heat5h_both_negative_high <- heat5h_both_negative[(abs(heat5h_both_negative$S_heat5h.log2FoldChange)>=4 & abs(heat5h_both_negative$D_heat5h.log2FoldChange)>=4),]
nrow(heat5h_both_negative_high)
#1
#heat5h_both_negative_high <- heat5h_both_negative[(abs(heat5h_both_negative$D_heat5h.log2FoldChange)>=4),]

#both positive
heat5h_both_positive <- heat5h_both_DGE[(heat5h_both_DGE$S_heat5h.log2FoldChange>0 & heat5h_both_DGE$D_heat5h.log2FoldChange>0),]
nrow(heat5h_both_positive)
#561
heat5h_both_positive_high <- heat5h_both_positive[(abs(heat5h_both_positive$S_heat5h.log2FoldChange)>=4 & abs(heat5h_both_positive$D_heat5h.log2FoldChange)>=4),]
nrow(heat5h_both_positive_high)
#0

#negative in S but positive in D
heat5h_pS_nD <- heat5h_both_DGE[(heat5h_both_DGE$S_heat5h.log2FoldChange>0 & heat5h_both_DGE$D_heat5h.log2FoldChange<0),]
nrow(heat5h_pS_nD)
#507
heat5h_pS_nD_high <- heat5h_pS_nD[(abs(heat5h_pS_nD$S_heat5h.log2FoldChange)>=4 & abs(heat5h_pS_nD$D_heat5h.log2FoldChange)>=4),]
nrow(heat5h_pS_nD_high)
#0
heat5h_nS_pD <- heat5h_both_DGE[(heat5h_both_DGE$S_heat5h.log2FoldChange<0 & heat5h_both_DGE$D_heat5h.log2FoldChange>0),]
nrow(heat5h_nS_pD)
#500

heat5h_nS_pD_high <- heat5h_nS_pD[(abs(heat5h_nS_pD$S_heat5h.log2FoldChange)>=4 & abs(heat5h_nS_pD$D_heat5h.log2FoldChange)>=4),]
nrow(heat5h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- heat5h_both_DGE$S_heat5h.log2FoldChange/heat5h_both_DGE$D_heat5h.log2FoldChange
heat5h_S_D_plasticity <- heat5h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(heat5h_S_D_plasticity)
#194 check the GO for those genes
heat5h_S_D_plasticity_high <- heat5h_S_D_plasticity[(abs(heat5h_S_D_plasticity$S_heat5h.log2FoldChange)>=4 & abs(heat5h_S_D_plasticity$D_heat5h.log2FoldChange)>=4),]
nrow(heat5h_S_D_plasticity_high)
#0

#make plot for the heat5h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat5h_DEG_S_D.pdf",width = 8,height = 8)
plot(heat5h_both_negative$S_heat5h.log2FoldChange,heat5h_both_negative$D_heat5h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-5h")
points(heat5h_both_positive$S_heat5h.log2FoldChange,heat5h_both_positive$D_heat5h.log2FoldChange,col="green")
points(heat5h_pS_nD$S_heat5h.log2FoldChange,heat5h_pS_nD$D_heat5h.log2FoldChange,col="blue")
points(heat5h_nS_pD$S_heat5h.log2FoldChange,heat5h_nS_pD$D_heat5h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat5h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(heat5h_both_negative_high$S_heat5h.log2FoldChange,heat5h_both_negative_high$D_heat5h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-5h")
points(heat5h_both_positive_high$S_heat5h.log2FoldChange,heat5h_both_positive_high$D_heat5h.log2FoldChange,col="green")
points(heat5h_pS_nD_high$S_heat5h.log2FoldChange,heat5h_pS_nD_high$D_heat5h.log2FoldChange,col="blue")
points(heat5h_nS_pD_high$S_heat5h.log2FoldChange,heat5h_nS_pD_high$D_heat5h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


heat5h_S_unique <- heat5h[(heat5h$S_heat5h.padj<0.05 & heat5h$D_heat5h.padj>0.05),]
nrow(heat5h_S_unique)
#3646
heat5h_S_unique_high <- heat5h_S_unique[(abs(heat5h_S_unique$S_heat5h.log2FoldChange)>=4),]
nrow(heat5h_S_unique_high)
#94
heat5h_D_unique <- heat5h[(heat5h$S_heat5h.padj>0.05 & heat5h$D_heat5h.padj<0.05),]
nrow(heat5h_D_unique)
#3208
heat5h_D_unique_high <- heat5h_D_unique[(abs(heat5h_D_unique$D_heat5h.log2FoldChange)>=4),]
nrow(heat5h_D_unique_high)
#106
head(heat5h_D_unique_high)

heat5h_both_NOT_DGE <- heat5h[(heat5h$S_heat5h.padj>0.05 & heat5h$D_heat5h.padj>0.05),]
nrow(heat5h_both_NOT_DGE)
#6142

#Heat10h
heat10h <- all[,c("Syl_ID","S_heat10h.log2FoldChange","S_heat10h.padj","distID","D_heat10h.log2FoldChange","D_heat10h.padj")]
head(heat10h)
nrow(heat10h)
#[1] 20735
heat10h <- na.omit(heat10h)
nrow(heat10h)
#15639
heat10h_both_DGE <- heat10h[(heat10h$S_heat10h.padj<0.05 & heat10h$D_heat10h.padj<0.05),]
nrow(heat10h_both_DGE)
#2401


#both negative
heat10h_both_negative <- heat10h_both_DGE[(heat10h_both_DGE$S_heat10h.log2FoldChange<0 & heat10h_both_DGE$D_heat10h.log2FoldChange<0),]
nrow(heat10h_both_negative)
#503
heat10h_both_negative_high <- heat10h_both_negative[(abs(heat10h_both_negative$S_heat10h.log2FoldChange)>=4 & abs(heat10h_both_negative$D_heat10h.log2FoldChange)>=4),]
nrow(heat10h_both_negative_high)
#1

#both positive
heat10h_both_positive <- heat10h_both_DGE[(heat10h_both_DGE$S_heat10h.log2FoldChange>0 & heat10h_both_DGE$D_heat10h.log2FoldChange>0),]
nrow(heat10h_both_positive)
#710
heat10h_both_positive_high <- heat10h_both_positive[(abs(heat10h_both_positive$S_heat10h.log2FoldChange)>=4 & abs(heat10h_both_positive$D_heat10h.log2FoldChange)>=4),]
nrow(heat10h_both_positive_high)
#1

#negative in S but positive in D
heat10h_pS_nD <- heat10h_both_DGE[(heat10h_both_DGE$S_heat10h.log2FoldChange>0 & heat10h_both_DGE$D_heat10h.log2FoldChange<0),]
nrow(heat10h_pS_nD)
#565
heat10h_pS_nD_high <- heat10h_pS_nD[(abs(heat10h_pS_nD$S_heat10h.log2FoldChange)>=4 & abs(heat10h_pS_nD$D_heat10h.log2FoldChange)>=4),]
nrow(heat10h_pS_nD_high)
#0
heat10h_nS_pD <- heat10h_both_DGE[(heat10h_both_DGE$S_heat10h.log2FoldChange<0 & heat10h_both_DGE$D_heat10h.log2FoldChange>0),]
nrow(heat10h_nS_pD)
#623

heat10h_nS_pD_high <- heat10h_nS_pD[(abs(heat10h_nS_pD$S_heat10h.log2FoldChange)>=4 & abs(heat10h_nS_pD$D_heat10h.log2FoldChange)>=4),]
nrow(heat10h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- heat10h_both_DGE$S_heat10h.log2FoldChange/heat10h_both_DGE$D_heat10h.log2FoldChange
heat10h_S_D_plasticity <- heat10h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(heat10h_S_D_plasticity)
#205 check the GO for those genes
heat10h_S_D_plasticity_high <- heat10h_S_D_plasticity[(abs(heat10h_S_D_plasticity$S_heat10h.log2FoldChange)>=4 & abs(heat10h_S_D_plasticity$D_heat10h.log2FoldChange)>=4),]
nrow(heat10h_S_D_plasticity_high)
#0

#make plot for the heat10h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat10h_DEG_S_D.pdf",width = 8,height = 8)
plot(heat10h_both_negative$S_heat10h.log2FoldChange,heat10h_both_negative$D_heat10h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-10h")
points(heat10h_both_positive$S_heat10h.log2FoldChange,heat10h_both_positive$D_heat10h.log2FoldChange,col="green")
points(heat10h_pS_nD$S_heat10h.log2FoldChange,heat10h_pS_nD$D_heat10h.log2FoldChange,col="blue")
points(heat10h_nS_pD$S_heat10h.log2FoldChange,heat10h_nS_pD$D_heat10h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat10h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(heat10h_both_negative_high$S_heat10h.log2FoldChange,heat10h_both_negative_high$D_heat10h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-10h")
points(heat10h_both_positive_high$S_heat10h.log2FoldChange,heat10h_both_positive_high$D_heat10h.log2FoldChange,col="green")
points(heat10h_pS_nD_high$S_heat10h.log2FoldChange,heat10h_pS_nD_high$D_heat10h.log2FoldChange,col="blue")
points(heat10h_nS_pD_high$S_heat10h.log2FoldChange,heat10h_nS_pD_high$D_heat10h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()

heat10h_S_unique <- heat10h[(heat10h$S_heat10h.padj<0.05 & heat10h$D_heat10h.padj>0.05),]
nrow(heat10h_S_unique)
#3950
heat10h_S_unique_high <- heat10h_S_unique[(abs(heat10h_S_unique$S_heat10h.log2FoldChange)>=4),]
nrow(heat10h_S_unique_high)
#140
heat10h_D_unique <- heat10h[(heat10h$S_heat10h.padj>0.05 & heat10h$D_heat10h.padj<0.05),]
nrow(heat10h_D_unique)
#3368
heat10h_D_unique_high <- heat10h_D_unique[(abs(heat10h_D_unique$D_heat10h.log2FoldChange)>=4),]
nrow(heat10h_D_unique_high)
#146
head(heat10h_D_unique_high)

heat10h_both_NOT_DGE <- heat10h[(heat10h$S_heat10h.padj>0.05 & heat10h$D_heat10h.padj>0.05),]
nrow(heat10h_both_NOT_DGE)
#5920


#24h:
#Heat24h
heat24h <- all[,c("Syl_ID","S_heat24h.log2FoldChange","S_heat24h.padj","distID","D_heat24h.log2FoldChange","D_heat24h.padj")]
head(heat24h)
nrow(heat24h)
#[1] 20735
heat24h <- na.omit(heat24h)
nrow(heat24h)
#16231
heat24h_both_DGE <- heat24h[(heat24h$S_heat24h.padj<0.05 & heat24h$D_heat24h.padj<0.05),]
nrow(heat24h_both_DGE)
#3586


#both negative
heat24h_both_negative <- heat24h_both_DGE[(heat24h_both_DGE$S_heat24h.log2FoldChange<0 & heat24h_both_DGE$D_heat24h.log2FoldChange<0),]
nrow(heat24h_both_negative)
#810
heat24h_both_negative_high <- heat24h_both_negative[(abs(heat24h_both_negative$S_heat24h.log2FoldChange)>=4 & abs(heat24h_both_negative$D_heat24h.log2FoldChange)>=4),]
nrow(heat24h_both_negative_high)
#0
#both positive
heat24h_both_positive <- heat24h_both_DGE[(heat24h_both_DGE$S_heat24h.log2FoldChange>0 & heat24h_both_DGE$D_heat24h.log2FoldChange>0),]
nrow(heat24h_both_positive)
#990
heat24h_both_positive_high <- heat24h_both_positive[(abs(heat24h_both_positive$S_heat24h.log2FoldChange)>=4 & abs(heat24h_both_positive$D_heat24h.log2FoldChange)>=4),]
nrow(heat24h_both_positive_high)
#4

#negative in S but positive in D
heat24h_pS_nD <- heat24h_both_DGE[(heat24h_both_DGE$S_heat24h.log2FoldChange>0 & heat24h_both_DGE$D_heat24h.log2FoldChange<0),]
nrow(heat24h_pS_nD)
#845
heat24h_pS_nD_high <- heat24h_pS_nD[(abs(heat24h_pS_nD$S_heat24h.log2FoldChange)>=4 & abs(heat24h_pS_nD$D_heat24h.log2FoldChange)>=4),]
nrow(heat24h_pS_nD_high)
#2
heat24h_nS_pD <- heat24h_both_DGE[(heat24h_both_DGE$S_heat24h.log2FoldChange<0 & heat24h_both_DGE$D_heat24h.log2FoldChange>0),]
nrow(heat24h_nS_pD)
#941

heat24h_nS_pD_high <- heat24h_nS_pD[(abs(heat24h_nS_pD$S_heat24h.log2FoldChange)>=4 & abs(heat24h_nS_pD$D_heat24h.log2FoldChange)>=4),]
nrow(heat24h_nS_pD_high)
#2

#look for the same direct and same fold changes
ratio <- heat24h_both_DGE$S_heat24h.log2FoldChange/heat24h_both_DGE$D_heat24h.log2FoldChange
heat24h_S_D_plasticity <- heat24h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(heat24h_S_D_plasticity)
#301 check the GO for those genes
heat24h_S_D_plasticity_high <- heat24h_S_D_plasticity[(abs(heat24h_S_D_plasticity$S_heat24h.log2FoldChange)>=4 & abs(heat24h_S_D_plasticity$D_heat24h.log2FoldChange)>=4),]
nrow(heat24h_S_D_plasticity_high)
#1

#make plot for the heat24h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat24h_DEG_S_D.pdf",width = 8,height = 8)
plot(heat24h_both_negative$S_heat24h.log2FoldChange,heat24h_both_negative$D_heat24h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-24h")
points(heat24h_both_positive$S_heat24h.log2FoldChange,heat24h_both_positive$D_heat24h.log2FoldChange,col="green")
points(heat24h_pS_nD$S_heat24h.log2FoldChange,heat24h_pS_nD$D_heat24h.log2FoldChange,col="blue")
points(heat24h_nS_pD$S_heat24h.log2FoldChange,heat24h_nS_pD$D_heat24h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat24h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(heat24h_both_negative_high$S_heat24h.log2FoldChange,heat24h_both_negative_high$D_heat24h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="Heat-24h")
points(heat24h_both_positive_high$S_heat24h.log2FoldChange,heat24h_both_positive_high$D_heat24h.log2FoldChange,col="green")
points(heat24h_pS_nD_high$S_heat24h.log2FoldChange,heat24h_pS_nD_high$D_heat24h.log2FoldChange,col="blue")
points(heat24h_nS_pD_high$S_heat24h.log2FoldChange,heat24h_nS_pD_high$D_heat24h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#3233


heat24h_S_unique <- heat24h[(heat24h$S_heat24h.padj<0.05 & heat24h$D_heat24h.padj>0.05),]
nrow(heat24h_S_unique)
#4374
heat24h_S_unique_high <- heat24h_S_unique[(abs(heat24h_S_unique$S_heat24h.log2FoldChange)>=4),]
nrow(heat24h_S_unique_high)
#152
heat24h_D_unique <- heat24h[(heat24h$S_heat24h.padj>0.05 & heat24h$D_heat24h.padj<0.05),]
nrow(heat24h_D_unique)
#3679
heat24h_D_unique_high <- heat24h_D_unique[(abs(heat24h_D_unique$D_heat24h.log2FoldChange)>=4),]
nrow(heat24h_D_unique_high)
#170

heat24h_both_NOT_DGE <- heat24h[(heat24h$S_heat24h.padj>0.05 & heat24h$D_heat24h.padj>0.05),]
nrow(heat24h_both_NOT_DGE)
#4529

#plot the barplot for the different categories for heat stress;

# Grouped Bar Plot
heat_plasticity <- c(nrow(heat1h_S_D_plasticity)/nrow(heat1h),
                     nrow(heat2h_S_D_plasticity)/nrow(heat2h),
                     nrow(heat5h_S_D_plasticity)/nrow(heat5h),
                     nrow(heat10h_S_D_plasticity)/nrow(heat10h),
                     nrow(heat24h_S_D_plasticity)/nrow(heat24h)
)

#high
heat_plasticity_high <- c(nrow(heat1h_S_D_plasticity_high)/nrow(heat1h),
                     nrow(heat2h_S_D_plasticity_high)/nrow(heat2h),
                     nrow(heat5h_S_D_plasticity_high)/nrow(heat5h),
                     nrow(heat10h_S_D_plasticity_high)/nrow(heat10h),
                     nrow(heat24h_S_D_plasticity_high)/nrow(heat24h)
)

heat_same_direction <- c(((nrow(heat1h_both_negative)+nrow(heat1h_both_positive)))/nrow(heat1h),
                         ((nrow(heat2h_both_negative)+nrow(heat2h_both_positive)))/nrow(heat2h),
                         ((nrow(heat5h_both_negative)+nrow(heat5h_both_positive)))/nrow(heat5h),
                         ((nrow(heat10h_both_negative)+nrow(heat10h_both_positive)))/nrow(heat10h),
                         ((nrow(heat24h_both_negative)+nrow(heat24h_both_positive)))/nrow(heat24h)
)

#high
heat_same_direction_high <- c(((nrow(heat1h_both_negative_high)+nrow(heat1h_both_positive_high)))/nrow(heat1h),
                         ((nrow(heat2h_both_negative_high)+nrow(heat2h_both_positive_high)))/nrow(heat2h),
                         ((nrow(heat5h_both_negative_high)+nrow(heat5h_both_positive_high)))/nrow(heat5h),
                         ((nrow(heat10h_both_negative_high)+nrow(heat10h_both_positive_high)))/nrow(heat10h),
                         ((nrow(heat24h_both_negative_high)+nrow(heat24h_both_positive_high)))/nrow(heat24h)
)

heat_opposite_direction <- c(((nrow(heat1h_pS_nD)+nrow(heat1h_nS_pD)))/nrow(heat1h),
                             ((nrow(heat2h_pS_nD)+nrow(heat2h_nS_pD)))/nrow(heat2h),
                             ((nrow(heat5h_pS_nD)+nrow(heat5h_nS_pD)))/nrow(heat5h),
                             ((nrow(heat10h_pS_nD)+nrow(heat10h_nS_pD)))/nrow(heat10h),
                             ((nrow(heat24h_pS_nD)+nrow(heat24h_nS_pD)))/nrow(heat24h)
)
#high
heat_opposite_direction_high <- c(((nrow(heat1h_pS_nD_high)+nrow(heat1h_nS_pD_high)))/nrow(heat1h),
                             ((nrow(heat2h_pS_nD_high)+nrow(heat2h_nS_pD_high)))/nrow(heat2h),
                             ((nrow(heat5h_pS_nD_high)+nrow(heat5h_nS_pD_high)))/nrow(heat5h),
                             ((nrow(heat10h_pS_nD_high)+nrow(heat10h_nS_pD_high)))/nrow(heat10h),
                             ((nrow(heat24h_pS_nD_high)+nrow(heat24h_nS_pD_high)))/nrow(heat24h)
)

heat_S_unique <- c(nrow(heat1h_S_unique)/nrow(heat1h),
                   nrow(heat2h_S_unique)/nrow(heat2h),
                   nrow(heat5h_S_unique)/nrow(heat5h),
                   nrow(heat10h_S_unique)/nrow(heat10h),
                   nrow(heat24h_S_unique)/nrow(heat24h)
)

#high
heat_S_unique_high <- c(nrow(heat1h_S_unique_high)/nrow(heat1h),
                   nrow(heat2h_S_unique_high)/nrow(heat2h),
                   nrow(heat5h_S_unique_high)/nrow(heat5h),
                   nrow(heat10h_S_unique_high)/nrow(heat10h),
                   nrow(heat24h_S_unique_high)/nrow(heat24h)
)
###save the unique gene list for heat treatment:
heat1h_S_unique_high_gene <- heat1h_S_unique_high$Syl_ID
heat2h_S_unique_high_gene <- heat2h_S_unique_high$Syl_ID
heat5h_S_unique_high_gene <- heat5h_S_unique_high$Syl_ID
heat10h_S_unique_high_gene <- heat10h_S_unique_high$Syl_ID
heat24h_S_unique_high_gene <- heat24h_S_unique_high$Syl_ID
write.table(
        heat1h_S_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat1h_S_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat2h_S_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat2h_S_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat5h_S_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat5h_S_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat10h_S_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat10h_S_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat24h_S_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat24h_S_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
###
heat_D_unique <- c(nrow(heat1h_D_unique)/nrow(heat1h),
                   nrow(heat2h_D_unique)/nrow(heat2h),
                   nrow(heat5h_D_unique)/nrow(heat5h),
                   nrow(heat10h_D_unique)/nrow(heat10h),
                   nrow(heat24h_D_unique)/nrow(heat24h)
)
#high
heat_D_unique_high <- c(nrow(heat1h_D_unique_high)/nrow(heat1h),
                   nrow(heat2h_D_unique_high)/nrow(heat2h),
                   nrow(heat5h_D_unique_high)/nrow(heat5h),
                   nrow(heat10h_D_unique_high)/nrow(heat10h),
                   nrow(heat24h_D_unique_high)/nrow(heat24h)
)

###save the D unique gene list for heat treatment:
heat1h_D_unique_high_gene <- heat1h_D_unique_high$distID
heat2h_D_unique_high_gene <- heat2h_D_unique_high$distID
heat5h_D_unique_high_gene <- heat5h_D_unique_high$distID
heat10h_D_unique_high_gene <- heat10h_D_unique_high$distID
heat24h_D_unique_high_gene <- heat24h_D_unique_high$distID
write.table(
        heat1h_D_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat1h_D_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat2h_D_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat2h_D_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat5h_D_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat5h_D_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat10h_D_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat10h_D_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")
write.table(
        heat24h_D_unique_high_gene,
        file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat24h_D_unique_high_gene.txt",
        quote=F,
        row.names=F,
        sep="\t")

heat24h_both_NOT_DGE <- c(nrow(heat1h_both_NOT_DGE)/nrow(heat1h),
                          nrow(heat2h_both_NOT_DGE)/nrow(heat2h),
                          nrow(heat5h_both_NOT_DGE)/nrow(heat5h),
                          nrow(heat10h_both_NOT_DGE)/nrow(heat10h),
                          nrow(heat24h_both_NOT_DGE)/nrow(heat24h)
)

heat_catNb <- cbind(heat_plasticity,
                    heat_same_direction,
                    heat_opposite_direction,
                    heat_S_unique,
                    heat_D_unique,
                    heat24h_both_NOT_DGE)
#high
heat_catNb_high <- cbind(heat_plasticity_high,
                    heat_same_direction_high,
                    heat_opposite_direction_high,
                    heat_S_unique_high,
                    heat_D_unique_high)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat_DEG_S_D.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(heat_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
barplot(heat_catNb,beside=T,col=(colfunc(5)),main="heat",ylim=c(0,0.8),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni","Both_exp_no_sig"))
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/heat_DEG_S_D_high.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(heat_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
barplot(heat_catNb_high,beside=T,col=(colfunc(5)),main="heat",ylim=c(0,0.03),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni"))
dev.off()

######################*****************************************################
#Drought
#
#####################*****************************************################
#drought1h
drought1h <- all[,c("Syl_ID","S_drought1h.log2FoldChange","S_drought1h.padj","distID","D_drought1h.log2FoldChange","D_drought1h.padj")]
head(drought1h)
nrow(drought1h)
#[1] 20735
drought1h <- na.omit(drought1h)
nrow(drought1h)
#17716
drought1h_both_DGE <- drought1h[(drought1h$S_drought1h.padj<0.05 & drought1h$D_drought1h.padj<0.05),]
nrow(drought1h_both_DGE)
#2343


#both negative
drought1h_both_negative <- drought1h_both_DGE[(drought1h_both_DGE$S_drought1h.log2FoldChange<0 & drought1h_both_DGE$D_drought1h.log2FoldChange<0),]
nrow(drought1h_both_negative)
#383
drought1h_both_negative_high <- drought1h_both_negative[(abs(drought1h_both_negative$S_drought1h.log2FoldChange)>=4 & abs(drought1h_both_negative$D_drought1h.log2FoldChange)>=4),]
nrow(drought1h_both_negative_high)
#1
#both positive
drought1h_both_positive <- drought1h_both_DGE[(drought1h_both_DGE$S_drought1h.log2FoldChange>0 & drought1h_both_DGE$D_drought1h.log2FoldChange>0),]
nrow(drought1h_both_positive)
#909
drought1h_both_positive_high <- drought1h_both_positive[(abs(drought1h_both_positive$S_drought1h.log2FoldChange)>=4 & abs(drought1h_both_positive$D_drought1h.log2FoldChange)>=4),]
nrow(drought1h_both_positive_high)
#104

#negative in S but positive in D
drought1h_pS_nD <- drought1h_both_DGE[(drought1h_both_DGE$S_drought1h.log2FoldChange>0 & drought1h_both_DGE$D_drought1h.log2FoldChange<0),]
nrow(drought1h_pS_nD)
#664
drought1h_pS_nD_high <- drought1h_pS_nD[(abs(drought1h_pS_nD$S_drought1h.log2FoldChange)>=4 & abs(drought1h_pS_nD$D_drought1h.log2FoldChange)>=4),]
nrow(drought1h_pS_nD_high)
#0
drought1h_nS_pD <- drought1h_both_DGE[(drought1h_both_DGE$S_drought1h.log2FoldChange<0 & drought1h_both_DGE$D_drought1h.log2FoldChange>0),]
nrow(drought1h_nS_pD)
#387

drought1h_nS_pD_high <- drought1h_nS_pD[(abs(drought1h_nS_pD$S_drought1h.log2FoldChange)>=4 & abs(drought1h_nS_pD$D_drought1h.log2FoldChange)>=4),]
nrow(drought1h_nS_pD_high)
#3

#look for the same direct and same fold changes
ratio <- drought1h_both_DGE$S_drought1h.log2FoldChange/drought1h_both_DGE$D_drought1h.log2FoldChange
drought1h_S_D_plasticity <- drought1h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(drought1h_S_D_plasticity)
#261 check the GO for those genes
drought1h_S_D_plasticity_high <- drought1h_S_D_plasticity[(abs(drought1h_S_D_plasticity$S_drought1h.log2FoldChange)>=4 & abs(drought1h_S_D_plasticity$D_drought1h.log2FoldChange)>=4),]
nrow(drought1h_S_D_plasticity_high)
#47

#make plot for the drought1h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought1h_DEG_S_D.pdf",width = 8,height = 8)
plot(drought1h_both_negative$S_drought1h.log2FoldChange,drought1h_both_negative$D_drought1h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-1h")
points(drought1h_both_positive$S_drought1h.log2FoldChange,drought1h_both_positive$D_drought1h.log2FoldChange,col="green")
points(drought1h_pS_nD$S_drought1h.log2FoldChange,drought1h_pS_nD$D_drought1h.log2FoldChange,col="blue")
points(drought1h_nS_pD$S_drought1h.log2FoldChange,drought1h_nS_pD$D_drought1h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought1h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(drought1h_both_negative_high$S_drought1h.log2FoldChange,drought1h_both_negative_high$D_drought1h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-1h")
points(drought1h_both_positive_high$S_drought1h.log2FoldChange,drought1h_both_positive_high$D_drought1h.log2FoldChange,col="green")
points(drought1h_pS_nD_high$S_drought1h.log2FoldChange,drought1h_pS_nD_high$D_drought1h.log2FoldChange,col="blue")
points(drought1h_nS_pD_high$S_drought1h.log2FoldChange,drought1h_nS_pD_high$D_drought1h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


drought1h_S_unique <- drought1h[(drought1h$S_drought1h.padj<0.05 & drought1h$D_drought1h.padj>0.05),]
nrow(drought1h_S_unique)
#2454
drought1h_S_unique_high <- drought1h_S_unique[(abs(drought1h_S_unique$S_drought1h.log2FoldChange)>=4),]
nrow(drought1h_S_unique_high)
#198
drought1h_D_unique <- drought1h[(drought1h$S_drought1h.padj>0.05 & drought1h$D_drought1h.padj<0.05),]
nrow(drought1h_D_unique)
#5110
drought1h_D_unique_high <- drought1h_D_unique[(abs(drought1h_D_unique$D_drought1h.log2FoldChange>=4)),]
nrow(drought1h_D_unique_high)
#196

drought1h_both_NOT_DGE <- drought1h[(drought1h$S_drought1h.padj>0.05 & drought1h$D_drought1h.padj>0.05),]
nrow(drought1h_both_NOT_DGE)
#7809

#drought2h
drought2h <- all[,c("Syl_ID","S_drought2h.log2FoldChange","S_drought2h.padj","distID","D_drought2h.log2FoldChange","D_drought2h.padj")]
head(drought2h)
nrow(drought2h)
#[1] 20735
drought2h <- na.omit(drought2h)
nrow(drought2h)
#18105
drought2h_both_DGE <- drought2h[(drought2h$S_drought2h.padj<0.05 & drought2h$D_drought2h.padj<0.05),]
nrow(drought2h_both_DGE)
#1642


#both negative
drought2h_both_negative <- drought2h_both_DGE[(drought2h_both_DGE$S_drought2h.log2FoldChange<0 & drought2h_both_DGE$D_drought2h.log2FoldChange<0),]
nrow(drought2h_both_negative)
#484
drought2h_both_negative_high <- drought2h_both_negative[(abs(drought2h_both_negative$S_drought2h.log2FoldChange)>=4 & abs(drought2h_both_negative$D_drought2h.log2FoldChange)>=4),]
nrow(drought2h_both_negative_high)
#6

#both positive
drought2h_both_positive <- drought2h_both_DGE[(drought2h_both_DGE$S_drought2h.log2FoldChange>0 & drought2h_both_DGE$D_drought2h.log2FoldChange>0),]
nrow(drought2h_both_positive)
#924
drought2h_both_positive_high <- drought2h_both_positive[(abs(drought2h_both_positive$S_drought2h.log2FoldChange)>=4 & abs(drought2h_both_positive$D_drought2h.log2FoldChange)>=4),]
nrow(drought2h_both_positive_high)
#120

#negative in S but positive in D
drought2h_pS_nD <- drought2h_both_DGE[(drought2h_both_DGE$S_drought2h.log2FoldChange>0 & drought2h_both_DGE$D_drought2h.log2FoldChange<0),]
nrow(drought2h_pS_nD)
#147
drought2h_pS_nD_high <- drought2h_pS_nD[(abs(drought2h_pS_nD$S_drought2h.log2FoldChange)>=4 & abs(drought2h_pS_nD$D_drought2h.log2FoldChange)>=4),]
nrow(drought2h_pS_nD_high)
#1
drought2h_nS_pD <- drought2h_both_DGE[(drought2h_both_DGE$S_drought2h.log2FoldChange<0 & drought2h_both_DGE$D_drought2h.log2FoldChange>0),]
nrow(drought2h_nS_pD)
#87

drought2h_nS_pD_high <- drought2h_nS_pD[(abs(drought2h_nS_pD$S_drought2h.log2FoldChange)>=4 & abs(drought2h_nS_pD$D_drought2h.log2FoldChange)>=4),]
nrow(drought2h_nS_pD_high)
#1

#look for the same direct and same fold changes
ratio <- drought2h_both_DGE$S_drought2h.log2FoldChange/drought2h_both_DGE$D_drought2h.log2FoldChange
drought2h_S_D_plasticity <- drought2h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(drought2h_S_D_plasticity)
#338 check the GO for those genes
drought2h_S_D_plasticity_high <- drought2h_S_D_plasticity[(abs(drought2h_S_D_plasticity$S_drought2h.log2FoldChange)>=4 & abs(drought2h_S_D_plasticity$D_drought2h.log2FoldChange)>=4),]
nrow(drought2h_S_D_plasticity_high)
#53

#make plot for the drought2h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought2h_DEG_S_D.pdf",width = 8,height = 8)
plot(drought2h_both_negative$S_drought2h.log2FoldChange,drought2h_both_negative$D_drought2h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-2h")
points(drought2h_both_positive$S_drought2h.log2FoldChange,drought2h_both_positive$D_drought2h.log2FoldChange,col="green")
points(drought2h_pS_nD$S_drought2h.log2FoldChange,drought2h_pS_nD$D_drought2h.log2FoldChange,col="blue")
points(drought2h_nS_pD$S_drought2h.log2FoldChange,drought2h_nS_pD$D_drought2h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought2h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(drought2h_both_negative_high$S_drought2h.log2FoldChange,drought2h_both_negative_high$D_drought2h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-2h")
points(drought2h_both_positive_high$S_drought2h.log2FoldChange,drought2h_both_positive_high$D_drought2h.log2FoldChange,col="green")
points(drought2h_pS_nD_high$S_drought2h.log2FoldChange,drought2h_pS_nD_high$D_drought2h.log2FoldChange,col="blue")
points(drought2h_nS_pD_high$S_drought2h.log2FoldChange,drought2h_nS_pD_high$D_drought2h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


drought2h_S_unique <- drought2h[(drought2h$S_drought2h.padj<0.05 & drought2h$D_drought2h.padj>0.05),]
nrow(drought2h_S_unique)
#6281
drought2h_S_unique_high <- drought2h_S_unique[(abs(drought2h_S_unique$S_drought2h.log2FoldChange)>=4),]
nrow(drought2h_S_unique_high)
#430
drought2h_D_unique <- drought2h[(drought2h$S_drought2h.padj>0.05 & drought2h$D_drought2h.padj<0.05),]
nrow(drought2h_D_unique)
#1031
drought2h_D_unique_high <- drought2h_D_unique[(abs(drought2h_D_unique$D_drought2h.log2FoldChange>=4)),]
nrow(drought2h_D_unique_high)
#93

drought2h_both_NOT_DGE <- drought2h[(drought2h$S_drought2h.padj>0.05 & drought2h$D_drought2h.padj>0.05),]
nrow(drought2h_both_NOT_DGE)
#9151

#5h:
drought5h <- all[,c("Syl_ID","S_drought5h.log2FoldChange","S_drought5h.padj","distID","D_drought5h.log2FoldChange","D_drought5h.padj")]
head(drought5h)
nrow(drought5h)
#[1] 20735
drought5h <- na.omit(drought5h)
nrow(drought5h)
#18190
drought5h_both_DGE <- drought5h[(drought5h$S_drought5h.padj<0.05 & drought5h$D_drought5h.padj<0.05),]
nrow(drought5h_both_DGE)
#5450


#both negative
drought5h_both_negative <- drought5h_both_DGE[(drought5h_both_DGE$S_drought5h.log2FoldChange<0 & drought5h_both_DGE$D_drought5h.log2FoldChange<0),]
nrow(drought5h_both_negative)
#2052
drought5h_both_negative_high <- drought5h_both_negative[(abs(drought5h_both_negative$S_drought5h.log2FoldChange)>=4 & abs(drought5h_both_negative$D_drought5h.log2FoldChange)>=4),]
nrow(drought5h_both_negative_high)
#43
#drought5h_both_negative_high <- drought5h_both_negative[(abs(drought5h_both_negative$D_drought5h.log2FoldChange)>=4),]

#both positive
drought5h_both_positive <- drought5h_both_DGE[(drought5h_both_DGE$S_drought5h.log2FoldChange>0 & drought5h_both_DGE$D_drought5h.log2FoldChange>0),]
nrow(drought5h_both_positive)
#2473
drought5h_both_positive_high <- drought5h_both_positive[(abs(drought5h_both_positive$S_drought5h.log2FoldChange)>=4 & abs(drought5h_both_positive$D_drought5h.log2FoldChange)>=4),]
nrow(drought5h_both_positive_high)
#204

#negative in S but positive in D
drought5h_pS_nD <- drought5h_both_DGE[(drought5h_both_DGE$S_drought5h.log2FoldChange>0 & drought5h_both_DGE$D_drought5h.log2FoldChange<0),]
nrow(drought5h_pS_nD)
#614
drought5h_pS_nD_high <- drought5h_pS_nD[(abs(drought5h_pS_nD$S_drought5h.log2FoldChange)>=4 & abs(drought5h_pS_nD$D_drought5h.log2FoldChange)>=4),]
nrow(drought5h_pS_nD_high)
#9
drought5h_nS_pD <- drought5h_both_DGE[(drought5h_both_DGE$S_drought5h.log2FoldChange<0 & drought5h_both_DGE$D_drought5h.log2FoldChange>0),]
nrow(drought5h_nS_pD)
#311

drought5h_nS_pD_high <- drought5h_nS_pD[(abs(drought5h_nS_pD$S_drought5h.log2FoldChange)>=4 & abs(drought5h_nS_pD$D_drought5h.log2FoldChange)>=4),]
nrow(drought5h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- drought5h_both_DGE$S_drought5h.log2FoldChange/drought5h_both_DGE$D_drought5h.log2FoldChange
drought5h_S_D_plasticity <- drought5h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(drought5h_S_D_plasticity)
#1121 check the GO for those genes
drought5h_S_D_plasticity_high <- drought5h_S_D_plasticity[(abs(drought5h_S_D_plasticity$S_drought5h.log2FoldChange)>=4 & abs(drought5h_S_D_plasticity$D_drought5h.log2FoldChange)>=4),]
nrow(drought5h_S_D_plasticity_high)
#121

#make plot for the drought5h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought5h_DEG_S_D.pdf",width = 8,height = 8)
plot(drought5h_both_negative$S_drought5h.log2FoldChange,drought5h_both_negative$D_drought5h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-5h")
points(drought5h_both_positive$S_drought5h.log2FoldChange,drought5h_both_positive$D_drought5h.log2FoldChange,col="green")
points(drought5h_pS_nD$S_drought5h.log2FoldChange,drought5h_pS_nD$D_drought5h.log2FoldChange,col="blue")
points(drought5h_nS_pD$S_drought5h.log2FoldChange,drought5h_nS_pD$D_drought5h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought5h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(drought5h_both_negative_high$S_drought5h.log2FoldChange,drought5h_both_negative_high$D_drought5h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-5h")
points(drought5h_both_positive_high$S_drought5h.log2FoldChange,drought5h_both_positive_high$D_drought5h.log2FoldChange,col="green")
points(drought5h_pS_nD_high$S_drought5h.log2FoldChange,drought5h_pS_nD_high$D_drought5h.log2FoldChange,col="blue")
points(drought5h_nS_pD_high$S_drought5h.log2FoldChange,drought5h_nS_pD_high$D_drought5h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


drought5h_S_unique <- drought5h[(drought5h$S_drought5h.padj<0.05 & drought5h$D_drought5h.padj>0.05),]
nrow(drought5h_S_unique)
#6361
drought5h_S_unique_high <- drought5h_S_unique[(abs(drought5h_S_unique$S_drought5h.log2FoldChange)>=4),]
nrow(drought5h_S_unique_high)
#415
drought5h_D_unique <- drought5h[(drought5h$S_drought5h.padj>0.05 & drought5h$D_drought5h.padj<0.05),]
nrow(drought5h_D_unique)
#2375
drought5h_D_unique_high <- drought5h_D_unique[(abs(drought5h_D_unique$D_drought5h.log2FoldChange>=4)),]
nrow(drought5h_D_unique_high)
#102

drought5h_both_NOT_DGE <- drought5h[(drought5h$S_drought5h.padj>0.05 & drought5h$D_drought5h.padj>0.05),]
nrow(drought5h_both_NOT_DGE)
#4004

#drought10h
drought10h <- all[,c("Syl_ID","S_drought10h.log2FoldChange","S_drought10h.padj","distID","D_drought10h.log2FoldChange","D_drought10h.padj")]
head(drought10h)
nrow(drought10h)
#[1] 20735
drought10h <- na.omit(drought10h)
nrow(drought10h)
#18294
drought10h_both_DGE <- drought10h[(drought10h$S_drought10h.padj<0.05 & drought10h$D_drought10h.padj<0.05),]
nrow(drought10h_both_DGE)
#6492


#both negative
drought10h_both_negative <- drought10h_both_DGE[(drought10h_both_DGE$S_drought10h.log2FoldChange<0 & drought10h_both_DGE$D_drought10h.log2FoldChange<0),]
nrow(drought10h_both_negative)
#2313
drought10h_both_negative_high <- drought10h_both_negative[(abs(drought10h_both_negative$S_drought10h.log2FoldChange)>=4 & abs(drought10h_both_negative$D_drought10h.log2FoldChange)>=4),]
nrow(drought10h_both_negative_high)
#35

#both positive
drought10h_both_positive <- drought10h_both_DGE[(drought10h_both_DGE$S_drought10h.log2FoldChange>0 & drought10h_both_DGE$D_drought10h.log2FoldChange>0),]
nrow(drought10h_both_positive)
#3061
drought10h_both_positive_high <- drought10h_both_positive[(abs(drought10h_both_positive$S_drought10h.log2FoldChange)>=4 & abs(drought10h_both_positive$D_drought10h.log2FoldChange)>=4),]
nrow(drought10h_both_positive_high)
#301

#negative in S but positive in D
drought10h_pS_nD <- drought10h_both_DGE[(drought10h_both_DGE$S_drought10h.log2FoldChange>0 & drought10h_both_DGE$D_drought10h.log2FoldChange<0),]
nrow(drought10h_pS_nD)
#683
drought10h_pS_nD_high <- drought10h_pS_nD[(abs(drought10h_pS_nD$S_drought10h.log2FoldChange)>=4 & abs(drought10h_pS_nD$D_drought10h.log2FoldChange)>=4),]
nrow(drought10h_pS_nD_high)
#10
drought10h_nS_pD <- drought10h_both_DGE[(drought10h_both_DGE$S_drought10h.log2FoldChange<0 & drought10h_both_DGE$D_drought10h.log2FoldChange>0),]
nrow(drought10h_nS_pD)
#435

drought10h_nS_pD_high <- drought10h_nS_pD[(abs(drought10h_nS_pD$S_drought10h.log2FoldChange)>=4 & abs(drought10h_nS_pD$D_drought10h.log2FoldChange)>=4),]
nrow(drought10h_nS_pD_high)
#3

#look for the same direct and same fold changes
ratio <- drought10h_both_DGE$S_drought10h.log2FoldChange/drought10h_both_DGE$D_drought10h.log2FoldChange
drought10h_S_D_plasticity <- drought10h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(drought10h_S_D_plasticity)
#1296 check the GO for those genes
drought10h_S_D_plasticity_high <- drought10h_S_D_plasticity[(abs(drought10h_S_D_plasticity$S_drought10h.log2FoldChange)>=4 & abs(drought10h_S_D_plasticity$D_drought10h.log2FoldChange)>=4),]
nrow(drought10h_S_D_plasticity_high)
#155

#make plot for the drought10h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought10h_DEG_S_D.pdf",width = 8,height = 8)
plot(drought10h_both_negative$S_drought10h.log2FoldChange,drought10h_both_negative$D_drought10h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-10h")
points(drought10h_both_positive$S_drought10h.log2FoldChange,drought10h_both_positive$D_drought10h.log2FoldChange,col="green")
points(drought10h_pS_nD$S_drought10h.log2FoldChange,drought10h_pS_nD$D_drought10h.log2FoldChange,col="blue")
points(drought10h_nS_pD$S_drought10h.log2FoldChange,drought10h_nS_pD$D_drought10h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought10h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(drought10h_both_negative_high$S_drought10h.log2FoldChange,drought10h_both_negative_high$D_drought10h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-10h")
points(drought10h_both_positive_high$S_drought10h.log2FoldChange,drought10h_both_positive_high$D_drought10h.log2FoldChange,col="green")
points(drought10h_pS_nD_high$S_drought10h.log2FoldChange,drought10h_pS_nD_high$D_drought10h.log2FoldChange,col="blue")
points(drought10h_nS_pD_high$S_drought10h.log2FoldChange,drought10h_nS_pD_high$D_drought10h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()

drought10h_S_unique <- drought10h[(drought10h$S_drought10h.padj<0.05 & drought10h$D_drought10h.padj>0.05),]
nrow(drought10h_S_unique)
#5488
drought10h_S_unique_high <- drought10h_S_unique[(abs(drought10h_S_unique$S_drought10h.log2FoldChange)>=4),]
nrow(drought10h_S_unique_high)
#339
drought10h_D_unique <- drought10h[(drought10h$S_drought10h.padj>0.05 & drought10h$D_drought10h.padj<0.05),]
nrow(drought10h_D_unique)
#2684
drought10h_D_unique_high <- drought10h_D_unique[(abs(drought10h_D_unique$D_drought10h.log2FoldChange>=4)),]
nrow(drought10h_D_unique_high)
#177

drought10h_both_NOT_DGE <- drought10h[(drought10h$S_drought10h.padj>0.05 & drought10h$D_drought10h.padj>0.05),]
nrow(drought10h_both_NOT_DGE)
#3630


#24h:
#drought24h
drought24h <- all[,c("Syl_ID","S_drought24h.log2FoldChange","S_drought24h.padj","distID","D_drought24h.log2FoldChange","D_drought24h.padj")]
head(drought24h)
nrow(drought24h)
#[1] 20735
drought24h <- na.omit(drought24h)
nrow(drought24h)
#18190
drought24h_both_DGE <- drought24h[(drought24h$S_drought24h.padj<0.05 & drought24h$D_drought24h.padj<0.05),]
nrow(drought24h_both_DGE)
#6090


#both negative
drought24h_both_negative <- drought24h_both_DGE[(drought24h_both_DGE$S_drought24h.log2FoldChange<0 & drought24h_both_DGE$D_drought24h.log2FoldChange<0),]
nrow(drought24h_both_negative)
#2042
drought24h_both_negative_high <- drought24h_both_negative[(abs(drought24h_both_negative$S_drought24h.log2FoldChange)>=4 & abs(drought24h_both_negative$D_drought24h.log2FoldChange)>=4),]
nrow(drought24h_both_negative_high)
#11
#both positive
drought24h_both_positive <- drought24h_both_DGE[(drought24h_both_DGE$S_drought24h.log2FoldChange>0 & drought24h_both_DGE$D_drought24h.log2FoldChange>0),]
nrow(drought24h_both_positive)
#3198
drought24h_both_positive_high <- drought24h_both_positive[(abs(drought24h_both_positive$S_drought24h.log2FoldChange)>=4 & abs(drought24h_both_positive$D_drought24h.log2FoldChange)>=4),]
nrow(drought24h_both_positive_high)
#403

#negative in S but positive in D
drought24h_pS_nD <- drought24h_both_DGE[(drought24h_both_DGE$S_drought24h.log2FoldChange>0 & drought24h_both_DGE$D_drought24h.log2FoldChange<0),]
nrow(drought24h_pS_nD)
#474
drought24h_pS_nD_high <- drought24h_pS_nD[(abs(drought24h_pS_nD$S_drought24h.log2FoldChange)>=4 & abs(drought24h_pS_nD$D_drought24h.log2FoldChange)>=4),]
nrow(drought24h_pS_nD_high)
#2
drought24h_nS_pD <- drought24h_both_DGE[(drought24h_both_DGE$S_drought24h.log2FoldChange<0 & drought24h_both_DGE$D_drought24h.log2FoldChange>0),]
nrow(drought24h_nS_pD)
#376

drought24h_nS_pD_high <- drought24h_nS_pD[(abs(drought24h_nS_pD$S_drought24h.log2FoldChange)>=4 & abs(drought24h_nS_pD$D_drought24h.log2FoldChange)>=4),]
nrow(drought24h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- drought24h_both_DGE$S_drought24h.log2FoldChange/drought24h_both_DGE$D_drought24h.log2FoldChange
drought24h_S_D_plasticity <- drought24h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(drought24h_S_D_plasticity)
#1265check the GO for those genes
drought24h_S_D_plasticity_high <- drought24h_S_D_plasticity[(abs(drought24h_S_D_plasticity$S_drought24h.log2FoldChange)>=4 & abs(drought24h_S_D_plasticity$D_drought24h.log2FoldChange)>=4),]
nrow(drought24h_S_D_plasticity_high)
#186

#make plot for the drought24h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought24h_DEG_S_D.pdf",width = 8,height = 8)
plot(drought24h_both_negative$S_drought24h.log2FoldChange,drought24h_both_negative$D_drought24h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-24h")
points(drought24h_both_positive$S_drought24h.log2FoldChange,drought24h_both_positive$D_drought24h.log2FoldChange,col="green")
points(drought24h_pS_nD$S_drought24h.log2FoldChange,drought24h_pS_nD$D_drought24h.log2FoldChange,col="blue")
points(drought24h_nS_pD$S_drought24h.log2FoldChange,drought24h_nS_pD$D_drought24h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought24h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(drought24h_both_negative_high$S_drought24h.log2FoldChange,drought24h_both_negative_high$D_drought24h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="drought-24h")
points(drought24h_both_positive_high$S_drought24h.log2FoldChange,drought24h_both_positive_high$D_drought24h.log2FoldChange,col="green")
points(drought24h_pS_nD_high$S_drought24h.log2FoldChange,drought24h_pS_nD_high$D_drought24h.log2FoldChange,col="blue")
points(drought24h_nS_pD_high$S_drought24h.log2FoldChange,drought24h_nS_pD_high$D_drought24h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()



drought24h_S_unique <- drought24h[(drought24h$S_drought24h.padj<0.05 & drought24h$D_drought24h.padj>0.05),]
nrow(drought24h_S_unique)
#6376
drought24h_S_unique_high <- drought24h_S_unique[(abs(drought24h_S_unique$S_drought24h.log2FoldChange)>=4),]
nrow(drought24h_S_unique_high)
#495
drought24h_D_unique <- drought24h[(drought24h$S_drought24h.padj>0.05 & drought24h$D_drought24h.padj<0.05),]
nrow(drought24h_D_unique)
#2017
drought24h_D_unique_high <- drought24h_D_unique[(abs(drought24h_D_unique$D_drought24h.log2FoldChange>=4)),]
nrow(drought24h_D_unique_high)
#158

drought24h_both_NOT_DGE <- drought24h[(drought24h$S_drought24h.padj>0.05 & drought24h$D_drought24h.padj>0.05),]
nrow(drought24h_both_NOT_DGE)
#3707

#plot the barplot for the different categories for drought stress;

# Grouped Bar Plot
drought_plasticity <- c(nrow(drought1h_S_D_plasticity)/nrow(drought1h),
                        nrow(drought2h_S_D_plasticity)/nrow(drought2h),
                        nrow(drought5h_S_D_plasticity)/nrow(drought5h),
                        nrow(drought10h_S_D_plasticity)/nrow(drought10h),
                        nrow(drought24h_S_D_plasticity)/nrow(drought24h)
)

#high
drought_plasticity_high <- c(nrow(drought1h_S_D_plasticity_high)/nrow(drought1h),
                             nrow(drought2h_S_D_plasticity_high)/nrow(drought2h),
                             nrow(drought5h_S_D_plasticity_high)/nrow(drought5h),
                             nrow(drought10h_S_D_plasticity_high)/nrow(drought10h),
                             nrow(drought24h_S_D_plasticity_high)/nrow(drought24h)
)

drought_same_direction <- c(((nrow(drought1h_both_negative)+nrow(drought1h_both_positive)))/nrow(drought1h),
                            ((nrow(drought2h_both_negative)+nrow(drought2h_both_positive)))/nrow(drought2h),
                            ((nrow(drought5h_both_negative)+nrow(drought5h_both_positive)))/nrow(drought5h),
                            ((nrow(drought10h_both_negative)+nrow(drought10h_both_positive)))/nrow(drought10h),
                            ((nrow(drought24h_both_negative)+nrow(drought24h_both_positive)))/nrow(drought24h)
)

#high
drought_same_direction_high <- c(((nrow(drought1h_both_negative_high)+nrow(drought1h_both_positive_high)))/nrow(drought1h),
                                 ((nrow(drought2h_both_negative_high)+nrow(drought2h_both_positive_high)))/nrow(drought2h),
                                 ((nrow(drought5h_both_negative_high)+nrow(drought5h_both_positive_high)))/nrow(drought5h),
                                 ((nrow(drought10h_both_negative_high)+nrow(drought10h_both_positive_high)))/nrow(drought10h),
                                 ((nrow(drought24h_both_negative_high)+nrow(drought24h_both_positive_high)))/nrow(drought24h)
)

drought_opposite_direction <- c(((nrow(drought1h_pS_nD)+nrow(drought1h_nS_pD)))/nrow(drought1h),
                                ((nrow(drought2h_pS_nD)+nrow(drought2h_nS_pD)))/nrow(drought2h),
                                ((nrow(drought5h_pS_nD)+nrow(drought5h_nS_pD)))/nrow(drought5h),
                                ((nrow(drought10h_pS_nD)+nrow(drought10h_nS_pD)))/nrow(drought10h),
                                ((nrow(drought24h_pS_nD)+nrow(drought24h_nS_pD)))/nrow(drought24h)
)
#high
drought_opposite_direction_high <- c(((nrow(drought1h_pS_nD_high)+nrow(drought1h_nS_pD_high)))/nrow(drought1h),
                                     ((nrow(drought2h_pS_nD_high)+nrow(drought2h_nS_pD_high)))/nrow(drought2h),
                                     ((nrow(drought5h_pS_nD_high)+nrow(drought5h_nS_pD_high)))/nrow(drought5h),
                                     ((nrow(drought10h_pS_nD_high)+nrow(drought10h_nS_pD_high)))/nrow(drought10h),
                                     ((nrow(drought24h_pS_nD_high)+nrow(drought24h_nS_pD_high)))/nrow(drought24h)
)

drought_S_unique <- c(nrow(drought1h_S_unique)/nrow(drought1h),
                      nrow(drought2h_S_unique)/nrow(drought2h),
                      nrow(drought5h_S_unique)/nrow(drought5h),
                      nrow(drought10h_S_unique)/nrow(drought10h),
                      nrow(drought24h_S_unique)/nrow(drought24h)
)

#high
drought_S_unique_high <- c(nrow(drought1h_S_unique_high)/nrow(drought1h),
                           nrow(drought2h_S_unique_high)/nrow(drought2h),
                           nrow(drought5h_S_unique_high)/nrow(drought5h),
                           nrow(drought10h_S_unique_high)/nrow(drought10h),
                           nrow(drought24h_S_unique_high)/nrow(drought24h)
)

drought_D_unique <- c(nrow(drought1h_D_unique)/nrow(drought1h),
                      nrow(drought2h_S_unique)/nrow(drought2h),
                      nrow(drought5h_S_unique)/nrow(drought5h),
                      nrow(drought10h_S_unique)/nrow(drought10h),
                      nrow(drought24h_S_unique)/nrow(drought24h)
)
#high
drought_D_unique_high <- c(nrow(drought1h_D_unique_high)/nrow(drought1h),
                           nrow(drought2h_S_unique_high)/nrow(drought2h),
                           nrow(drought5h_S_unique_high)/nrow(drought5h),
                           nrow(drought10h_S_unique_high)/nrow(drought10h),
                           nrow(drought24h_S_unique_high)/nrow(drought24h)
)

drought24h_both_NOT_DGE <- c(nrow(drought1h_both_NOT_DGE)/nrow(drought1h),
                             nrow(drought2h_both_NOT_DGE)/nrow(drought2h),
                             nrow(drought5h_both_NOT_DGE)/nrow(drought5h),
                             nrow(drought10h_both_NOT_DGE)/nrow(drought10h),
                             nrow(drought24h_both_NOT_DGE)/nrow(drought24h)
)

drought_catNb <- cbind(drought_plasticity,
                       drought_same_direction,
                       drought_opposite_direction,
                       drought_S_unique,
                       drought_D_unique,
                       drought24h_both_NOT_DGE)
#high
drought_catNb_high <- cbind(drought_plasticity_high,
                            drought_same_direction_high,
                            drought_opposite_direction_high,
                            drought_S_unique_high,
                            drought_D_unique_high)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought_DEG_S_D.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(drought_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
barplot(drought_catNb,beside=T,col=(colfunc(5)),main="drought",ylim=c(0,0.8),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni","Both_exp_no_sig"))
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/drought_DEG_S_D_high.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(drought_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
barplot(drought_catNb_high,beside=T,col=(colfunc(5)),main="drought",ylim=c(0,0.03),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni"))
dev.off()

######################*****************************************################
#salt
#
#####################*****************************************################
#salt1h
salt1h <- all[,c("Syl_ID","S_salt1h.log2FoldChange","S_salt1h.padj","distID","D_salt1h.log2FoldChange","D_salt1h.padj")]
head(salt1h)
nrow(salt1h)
#[1] 20735
salt1h <- na.omit(salt1h)
nrow(salt1h)
#16982
salt1h_both_DGE <- salt1h[(salt1h$S_salt1h.padj<0.05 & salt1h$D_salt1h.padj<0.05),]
nrow(salt1h_both_DGE)
#840


#both negative
salt1h_both_negative <- salt1h_both_DGE[(salt1h_both_DGE$S_salt1h.log2FoldChange<0 & salt1h_both_DGE$D_salt1h.log2FoldChange<0),]
nrow(salt1h_both_negative)
#67
salt1h_both_negative_high <- salt1h_both_negative[(abs(salt1h_both_negative$S_salt1h.log2FoldChange)>=4 & abs(salt1h_both_negative$D_salt1h.log2FoldChange)>=4),]
nrow(salt1h_both_negative_high)
#0
#both positive
salt1h_both_positive <- salt1h_both_DGE[(salt1h_both_DGE$S_salt1h.log2FoldChange>0 & salt1h_both_DGE$D_salt1h.log2FoldChange>0),]
nrow(salt1h_both_positive)
#136
salt1h_both_positive_high <- salt1h_both_positive[(abs(salt1h_both_positive$S_salt1h.log2FoldChange)>=4 & abs(salt1h_both_positive$D_salt1h.log2FoldChange)>=4),]
nrow(salt1h_both_positive_high)
#2

#negative in S but positive in D
salt1h_pS_nD <- salt1h_both_DGE[(salt1h_both_DGE$S_salt1h.log2FoldChange>0 & salt1h_both_DGE$D_salt1h.log2FoldChange<0),]
nrow(salt1h_pS_nD)
#557
salt1h_pS_nD_high <- salt1h_pS_nD[(abs(salt1h_pS_nD$S_salt1h.log2FoldChange)>=4 & abs(salt1h_pS_nD$D_salt1h.log2FoldChange)>=4),]
nrow(salt1h_pS_nD_high)
#8
salt1h_nS_pD <- salt1h_both_DGE[(salt1h_both_DGE$S_salt1h.log2FoldChange<0 & salt1h_both_DGE$D_salt1h.log2FoldChange>0),]
nrow(salt1h_nS_pD)
#80

salt1h_nS_pD_high <- salt1h_nS_pD[(abs(salt1h_nS_pD$S_salt1h.log2FoldChange)>=4 & abs(salt1h_nS_pD$D_salt1h.log2FoldChange)>=4),]
nrow(salt1h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- salt1h_both_DGE$S_salt1h.log2FoldChange/salt1h_both_DGE$D_salt1h.log2FoldChange
salt1h_S_D_plasticity <- salt1h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(salt1h_S_D_plasticity)
#30 check the GO for those genes
salt1h_S_D_plasticity_high <- salt1h_S_D_plasticity[(abs(salt1h_S_D_plasticity$S_salt1h.log2FoldChange)>=4 & abs(salt1h_S_D_plasticity$D_salt1h.log2FoldChange)>=4),]
nrow(salt1h_S_D_plasticity_high)
#1

#make plot for the salt1h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt1h_DEG_S_D.pdf",width = 8,height = 8)
plot(salt1h_both_negative$S_salt1h.log2FoldChange,salt1h_both_negative$D_salt1h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-1h")
points(salt1h_both_positive$S_salt1h.log2FoldChange,salt1h_both_positive$D_salt1h.log2FoldChange,col="green")
points(salt1h_pS_nD$S_salt1h.log2FoldChange,salt1h_pS_nD$D_salt1h.log2FoldChange,col="blue")
points(salt1h_nS_pD$S_salt1h.log2FoldChange,salt1h_nS_pD$D_salt1h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt1h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(salt1h_both_negative_high$S_salt1h.log2FoldChange,salt1h_both_negative_high$D_salt1h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-1h")
points(salt1h_both_positive_high$S_salt1h.log2FoldChange,salt1h_both_positive_high$D_salt1h.log2FoldChange,col="green")
points(salt1h_pS_nD_high$S_salt1h.log2FoldChange,salt1h_pS_nD_high$D_salt1h.log2FoldChange,col="blue")
points(salt1h_nS_pD_high$S_salt1h.log2FoldChange,salt1h_nS_pD_high$D_salt1h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


salt1h_S_unique <- salt1h[(salt1h$S_salt1h.padj<0.05 & salt1h$D_salt1h.padj>0.05),]
nrow(salt1h_S_unique)
#1668
salt1h_S_unique_high <- salt1h_S_unique[(abs(salt1h_S_unique$S_salt1h.log2FoldChange)>=4),]
nrow(salt1h_S_unique_high)
#104
salt1h_D_unique <- salt1h[(salt1h$S_salt1h.padj>0.05 & salt1h$D_salt1h.padj<0.05),]
nrow(salt1h_D_unique)
#3365
salt1h_D_unique_high <- salt1h_D_unique[(abs(salt1h_D_unique$D_salt1h.log2FoldChange>=4)),]
nrow(salt1h_D_unique_high)
#42

salt1h_both_NOT_DGE <- salt1h[(salt1h$S_salt1h.padj>0.05 & salt1h$D_salt1h.padj>0.05),]
nrow(salt1h_both_NOT_DGE)
#11109

#salt2h
salt2h <- all[,c("Syl_ID","S_salt2h.log2FoldChange","S_salt2h.padj","distID","D_salt2h.log2FoldChange","D_salt2h.padj")]
head(salt2h)
nrow(salt2h)
#[1] 20735
salt2h <- na.omit(salt2h)
nrow(salt2h)
#16480
salt2h_both_DGE <- salt2h[(salt2h$S_salt2h.padj<0.05 & salt2h$D_salt2h.padj<0.05),]
nrow(salt2h_both_DGE)
#336


#both negative
salt2h_both_negative <- salt2h_both_DGE[(salt2h_both_DGE$S_salt2h.log2FoldChange<0 & salt2h_both_DGE$D_salt2h.log2FoldChange<0),]
nrow(salt2h_both_negative)
#91
salt2h_both_negative_high <- salt2h_both_negative[(abs(salt2h_both_negative$S_salt2h.log2FoldChange)>=4 & abs(salt2h_both_negative$D_salt2h.log2FoldChange)>=4),]
nrow(salt2h_both_negative_high)
#0
#salt2h_both_negative_high <- salt2h_both_negative[(abs(salt2h_both_negative$D_salt2h.log2FoldChange)>=4),]

#both positive
salt2h_both_positive <- salt2h_both_DGE[(salt2h_both_DGE$S_salt2h.log2FoldChange>0 & salt2h_both_DGE$D_salt2h.log2FoldChange>0),]
nrow(salt2h_both_positive)
#43
salt2h_both_positive_high <- salt2h_both_positive[(abs(salt2h_both_positive$S_salt2h.log2FoldChange)>=4 & abs(salt2h_both_positive$D_salt2h.log2FoldChange)>=4),]
nrow(salt2h_both_positive_high)
#0

#negative in S but positive in D
salt2h_pS_nD <- salt2h_both_DGE[(salt2h_both_DGE$S_salt2h.log2FoldChange>0 & salt2h_both_DGE$D_salt2h.log2FoldChange<0),]
nrow(salt2h_pS_nD)
#163
salt2h_pS_nD_high <- salt2h_pS_nD[(abs(salt2h_pS_nD$S_salt2h.log2FoldChange)>=4 & abs(salt2h_pS_nD$D_salt2h.log2FoldChange)>=4),]
nrow(salt2h_pS_nD_high)
#1
salt2h_nS_pD <- salt2h_both_DGE[(salt2h_both_DGE$S_salt2h.log2FoldChange<0 & salt2h_both_DGE$D_salt2h.log2FoldChange>0),]
nrow(salt2h_nS_pD)
#39

salt2h_nS_pD_high <- salt2h_nS_pD[(abs(salt2h_nS_pD$S_salt2h.log2FoldChange)>=4 & abs(salt2h_nS_pD$D_salt2h.log2FoldChange)>=4),]
nrow(salt2h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- salt2h_both_DGE$S_salt2h.log2FoldChange/salt2h_both_DGE$D_salt2h.log2FoldChange
salt2h_S_D_plasticity <- salt2h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(salt2h_S_D_plasticity)
#25check the GO for those genes
salt2h_S_D_plasticity_high <- salt2h_S_D_plasticity[(abs(salt2h_S_D_plasticity$S_salt2h.log2FoldChange)>=4 & abs(salt2h_S_D_plasticity$D_salt2h.log2FoldChange)>=4),]
nrow(salt2h_S_D_plasticity_high)
#0

#make plot for the salt2h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt2h_DEG_S_D.pdf",width = 8,height = 8)
plot(salt2h_both_negative$S_salt2h.log2FoldChange,salt2h_both_negative$D_salt2h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-2h")
points(salt2h_both_positive$S_salt2h.log2FoldChange,salt2h_both_positive$D_salt2h.log2FoldChange,col="green")
points(salt2h_pS_nD$S_salt2h.log2FoldChange,salt2h_pS_nD$D_salt2h.log2FoldChange,col="blue")
points(salt2h_nS_pD$S_salt2h.log2FoldChange,salt2h_nS_pD$D_salt2h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt2h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(salt2h_both_negative_high$S_salt2h.log2FoldChange,salt2h_both_negative_high$D_salt2h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-2h")
points(salt2h_both_positive_high$S_salt2h.log2FoldChange,salt2h_both_positive_high$D_salt2h.log2FoldChange,col="green")
points(salt2h_pS_nD_high$S_salt2h.log2FoldChange,salt2h_pS_nD_high$D_salt2h.log2FoldChange,col="blue")
points(salt2h_nS_pD_high$S_salt2h.log2FoldChange,salt2h_nS_pD_high$D_salt2h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#3233


salt2h_S_unique <- salt2h[(salt2h$S_salt2h.padj<0.05 & salt2h$D_salt2h.padj>0.05),]
nrow(salt2h_S_unique)
#1092
salt2h_S_unique_high <- salt2h_S_unique[(abs(salt2h_S_unique$S_salt2h.log2FoldChange)>=4),]
nrow(salt2h_S_unique_high)
#36
salt2h_D_unique <- salt2h[(salt2h$S_salt2h.padj>0.05 & salt2h$D_salt2h.padj<0.05),]
nrow(salt2h_D_unique)
#3557
salt2h_D_unique_high <- salt2h_D_unique[(abs(salt2h_D_unique$D_salt2h.log2FoldChange>=4)),]
nrow(salt2h_D_unique_high)
#19

salt2h_both_NOT_DGE <- salt2h[(salt2h$S_salt2h.padj>0.05 & salt2h$D_salt2h.padj>0.05),]
nrow(salt2h_both_NOT_DGE)
#11495

#5h:
salt5h <- all[,c("Syl_ID","S_salt5h.log2FoldChange","S_salt5h.padj","distID","D_salt5h.log2FoldChange","D_salt5h.padj")]
head(salt5h)
nrow(salt5h)
#[1] 20735
salt5h <- na.omit(salt5h)
nrow(salt5h)
#15301
salt5h_both_DGE <- salt5h[(salt5h$S_salt5h.padj<0.05 & salt5h$D_salt5h.padj<0.05),]
nrow(salt5h_both_DGE)
#228


#both negative
salt5h_both_negative <- salt5h_both_DGE[(salt5h_both_DGE$S_salt5h.log2FoldChange<0 & salt5h_both_DGE$D_salt5h.log2FoldChange<0),]
nrow(salt5h_both_negative)
#80
salt5h_both_negative_high <- salt5h_both_negative[(abs(salt5h_both_negative$S_salt5h.log2FoldChange)>=4 & abs(salt5h_both_negative$D_salt5h.log2FoldChange)>=4),]
nrow(salt5h_both_negative_high)
#0
#salt5h_both_negative_high <- salt5h_both_negative[(abs(salt5h_both_negative$D_salt5h.log2FoldChange)>=4),]

#both positive
salt5h_both_positive <- salt5h_both_DGE[(salt5h_both_DGE$S_salt5h.log2FoldChange>0 & salt5h_both_DGE$D_salt5h.log2FoldChange>0),]
nrow(salt5h_both_positive)
#116
salt5h_both_positive_high <- salt5h_both_positive[(abs(salt5h_both_positive$S_salt5h.log2FoldChange)>=4 & abs(salt5h_both_positive$D_salt5h.log2FoldChange)>=4),]
nrow(salt5h_both_positive_high)
#1

#negative in S but positive in D
salt5h_pS_nD <- salt5h_both_DGE[(salt5h_both_DGE$S_salt5h.log2FoldChange>0 & salt5h_both_DGE$D_salt5h.log2FoldChange<0),]
nrow(salt5h_pS_nD)
#26
salt5h_pS_nD_high <- salt5h_pS_nD[(abs(salt5h_pS_nD$S_salt5h.log2FoldChange)>=4 & abs(salt5h_pS_nD$D_salt5h.log2FoldChange)>=4),]
nrow(salt5h_pS_nD_high)
#0
salt5h_nS_pD <- salt5h_both_DGE[(salt5h_both_DGE$S_salt5h.log2FoldChange<0 & salt5h_both_DGE$D_salt5h.log2FoldChange>0),]
nrow(salt5h_nS_pD)
#6

salt5h_nS_pD_high <- salt5h_nS_pD[(abs(salt5h_nS_pD$S_salt5h.log2FoldChange)>=4 & abs(salt5h_nS_pD$D_salt5h.log2FoldChange)>=4),]
nrow(salt5h_nS_pD_high)
#1

#look for the same direct and same fold changes
ratio <- salt5h_both_DGE$S_salt5h.log2FoldChange/salt5h_both_DGE$D_salt5h.log2FoldChange
salt5h_S_D_plasticity <- salt5h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(salt5h_S_D_plasticity)
#28 check the GO for those genes
salt5h_S_D_plasticity_high <- salt5h_S_D_plasticity[(abs(salt5h_S_D_plasticity$S_salt5h.log2FoldChange)>=4 & abs(salt5h_S_D_plasticity$D_salt5h.log2FoldChange)>=4),]
nrow(salt5h_S_D_plasticity_high)
#1

#make plot for the salt5h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt5h_DEG_S_D.pdf",width = 8,height = 8)
plot(salt5h_both_negative$S_salt5h.log2FoldChange,salt5h_both_negative$D_salt5h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-5h")
points(salt5h_both_positive$S_salt5h.log2FoldChange,salt5h_both_positive$D_salt5h.log2FoldChange,col="green")
points(salt5h_pS_nD$S_salt5h.log2FoldChange,salt5h_pS_nD$D_salt5h.log2FoldChange,col="blue")
points(salt5h_nS_pD$S_salt5h.log2FoldChange,salt5h_nS_pD$D_salt5h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt5h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(salt5h_both_negative_high$S_salt5h.log2FoldChange,salt5h_both_negative_high$D_salt5h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-5h")
points(salt5h_both_positive_high$S_salt5h.log2FoldChange,salt5h_both_positive_high$D_salt5h.log2FoldChange,col="green")
points(salt5h_pS_nD_high$S_salt5h.log2FoldChange,salt5h_pS_nD_high$D_salt5h.log2FoldChange,col="blue")
points(salt5h_nS_pD_high$S_salt5h.log2FoldChange,salt5h_nS_pD_high$D_salt5h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()


salt5h_S_unique <- salt5h[(salt5h$S_salt5h.padj<0.05 & salt5h$D_salt5h.padj>0.05),]
nrow(salt5h_S_unique)
#681
salt5h_S_unique_high <- salt5h_S_unique[(abs(salt5h_S_unique$S_salt5h.log2FoldChange)>=4),]
nrow(salt5h_S_unique_high)
#5
salt5h_D_unique <- salt5h[(salt5h$S_salt5h.padj>0.05 & salt5h$D_salt5h.padj<0.05),]
nrow(salt5h_D_unique)
#1876
salt5h_D_unique_high <- salt5h_D_unique[(abs(salt5h_D_unique$D_salt5h.log2FoldChange>=4)),]
nrow(salt5h_D_unique_high)
#89

salt5h_both_NOT_DGE <- salt5h[(salt5h$S_salt5h.padj>0.05 & salt5h$D_salt5h.padj>0.05),]
nrow(salt5h_both_NOT_DGE)
#12516

#salt10h
salt10h <- all[,c("Syl_ID","S_salt10h.log2FoldChange","S_salt10h.padj","distID","D_salt10h.log2FoldChange","D_salt10h.padj")]
head(salt10h)
nrow(salt10h)
#[1] 20735
salt10h <- na.omit(salt10h)
nrow(salt10h)
#16776
salt10h_both_DGE <- salt10h[(salt10h$S_salt10h.padj<0.05 & salt10h$D_salt10h.padj<0.05),]
nrow(salt10h_both_DGE)
#519


#both negative
salt10h_both_negative <- salt10h_both_DGE[(salt10h_both_DGE$S_salt10h.log2FoldChange<0 & salt10h_both_DGE$D_salt10h.log2FoldChange<0),]
nrow(salt10h_both_negative)
#121
salt10h_both_negative_high <- salt10h_both_negative[(abs(salt10h_both_negative$S_salt10h.log2FoldChange)>=4 & abs(salt10h_both_negative$D_salt10h.log2FoldChange)>=4),]
nrow(salt10h_both_negative_high)
#0

#both positive
salt10h_both_positive <- salt10h_both_DGE[(salt10h_both_DGE$S_salt10h.log2FoldChange>0 & salt10h_both_DGE$D_salt10h.log2FoldChange>0),]
nrow(salt10h_both_positive)
#302
salt10h_both_positive_high <- salt10h_both_positive[(abs(salt10h_both_positive$S_salt10h.log2FoldChange)>=4 & abs(salt10h_both_positive$D_salt10h.log2FoldChange)>=4),]
nrow(salt10h_both_positive_high)
#9

#negative in S but positive in D
salt10h_pS_nD <- salt10h_both_DGE[(salt10h_both_DGE$S_salt10h.log2FoldChange>0 & salt10h_both_DGE$D_salt10h.log2FoldChange<0),]
nrow(salt10h_pS_nD)
#81
salt10h_pS_nD_high <- salt10h_pS_nD[(abs(salt10h_pS_nD$S_salt10h.log2FoldChange)>=4 & abs(salt10h_pS_nD$D_salt10h.log2FoldChange)>=4),]
nrow(salt10h_pS_nD_high)
#0
salt10h_nS_pD <- salt10h_both_DGE[(salt10h_both_DGE$S_salt10h.log2FoldChange<0 & salt10h_both_DGE$D_salt10h.log2FoldChange>0),]
nrow(salt10h_nS_pD)
#15

salt10h_nS_pD_high <- salt10h_nS_pD[(abs(salt10h_nS_pD$S_salt10h.log2FoldChange)>=4 & abs(salt10h_nS_pD$D_salt10h.log2FoldChange)>=4),]
nrow(salt10h_nS_pD_high)
#0

#look for the same direct and same fold changes
ratio <- salt10h_both_DGE$S_salt10h.log2FoldChange/salt10h_both_DGE$D_salt10h.log2FoldChange
salt10h_S_D_plasticity <- salt10h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(salt10h_S_D_plasticity)
#57 check the GO for those genes
salt10h_S_D_plasticity_high <- salt10h_S_D_plasticity[(abs(salt10h_S_D_plasticity$S_salt10h.log2FoldChange)>=4 & abs(salt10h_S_D_plasticity$D_salt10h.log2FoldChange)>=4),]
nrow(salt10h_S_D_plasticity_high)
#4

#make plot for the salt10h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt10h_DEG_S_D.pdf",width = 8,height = 8)
plot(salt10h_both_negative$S_salt10h.log2FoldChange,salt10h_both_negative$D_salt10h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-10h")
points(salt10h_both_positive$S_salt10h.log2FoldChange,salt10h_both_positive$D_salt10h.log2FoldChange,col="green")
points(salt10h_pS_nD$S_salt10h.log2FoldChange,salt10h_pS_nD$D_salt10h.log2FoldChange,col="blue")
points(salt10h_nS_pD$S_salt10h.log2FoldChange,salt10h_nS_pD$D_salt10h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt10h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(salt10h_both_negative_high$S_salt10h.log2FoldChange,salt10h_both_negative_high$D_salt10h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-10h")
points(salt10h_both_positive_high$S_salt10h.log2FoldChange,salt10h_both_positive_high$D_salt10h.log2FoldChange,col="green")
points(salt10h_pS_nD_high$S_salt10h.log2FoldChange,salt10h_pS_nD_high$D_salt10h.log2FoldChange,col="blue")
points(salt10h_nS_pD_high$S_salt10h.log2FoldChange,salt10h_nS_pD_high$D_salt10h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()

salt10h_S_unique <- salt10h[(salt10h$S_salt10h.padj<0.05 & salt10h$D_salt10h.padj>0.05),]
nrow(salt10h_S_unique)
#1324
salt10h_S_unique_high <- salt10h_S_unique[(abs(salt10h_S_unique$S_salt10h.log2FoldChange)>=4),]
nrow(salt10h_S_unique_high)
#20
salt10h_D_unique <- salt10h[(salt10h$S_salt10h.padj>0.05 & salt10h$D_salt10h.padj<0.05),]
nrow(salt10h_D_unique)
#3349
salt10h_D_unique_high <- salt10h_D_unique[(abs(salt10h_D_unique$D_salt10h.log2FoldChange>=4)),]
nrow(salt10h_D_unique_high)
#155

salt10h_both_NOT_DGE <- salt10h[(salt10h$S_salt10h.padj>0.05 & salt10h$D_salt10h.padj>0.05),]
nrow(salt10h_both_NOT_DGE)
#11584


#24h:
#salt24h
salt24h <- all[,c("Syl_ID","S_salt24h.log2FoldChange","S_salt24h.padj","distID","D_salt24h.log2FoldChange","D_salt24h.padj")]
head(salt24h)
nrow(salt24h)
#[1] 20735
salt24h <- na.omit(salt24h)
nrow(salt24h)
#16776
salt24h_both_DGE <- salt24h[(salt24h$S_salt24h.padj<0.05 & salt24h$D_salt24h.padj<0.05),]
nrow(salt24h_both_DGE)
#1475


#both negative
salt24h_both_negative <- salt24h_both_DGE[(salt24h_both_DGE$S_salt24h.log2FoldChange<0 & salt24h_both_DGE$D_salt24h.log2FoldChange<0),]
nrow(salt24h_both_negative)
#721
salt24h_both_negative_high <- salt24h_both_negative[(abs(salt24h_both_negative$S_salt24h.log2FoldChange)>=4 & abs(salt24h_both_negative$D_salt24h.log2FoldChange)>=4),]
nrow(salt24h_both_negative_high)
#2
#both positive
salt24h_both_positive <- salt24h_both_DGE[(salt24h_both_DGE$S_salt24h.log2FoldChange>0 & salt24h_both_DGE$D_salt24h.log2FoldChange>0),]
nrow(salt24h_both_positive)
#441
salt24h_both_positive_high <- salt24h_both_positive[(abs(salt24h_both_positive$S_salt24h.log2FoldChange)>=4 & abs(salt24h_both_positive$D_salt24h.log2FoldChange)>=4),]
nrow(salt24h_both_positive_high)
#17

#negative in S but positive in D
salt24h_pS_nD <- salt24h_both_DGE[(salt24h_both_DGE$S_salt24h.log2FoldChange>0 & salt24h_both_DGE$D_salt24h.log2FoldChange<0),]
nrow(salt24h_pS_nD)
#143
salt24h_pS_nD_high <- salt24h_pS_nD[(abs(salt24h_pS_nD$S_salt24h.log2FoldChange)>=4 & abs(salt24h_pS_nD$D_salt24h.log2FoldChange)>=4),]
nrow(salt24h_pS_nD_high)
#0
salt24h_nS_pD <- salt24h_both_DGE[(salt24h_both_DGE$S_salt24h.log2FoldChange<0 & salt24h_both_DGE$D_salt24h.log2FoldChange>0),]
nrow(salt24h_nS_pD)
#170

salt24h_nS_pD_high <- salt24h_nS_pD[(abs(salt24h_nS_pD$S_salt24h.log2FoldChange)>=4 & abs(salt24h_nS_pD$D_salt24h.log2FoldChange)>=4),]
nrow(salt24h_nS_pD_high)
#2

#look for the same direct and same fold changes
ratio <- salt24h_both_DGE$S_salt24h.log2FoldChange/salt24h_both_DGE$D_salt24h.log2FoldChange
salt24h_S_D_plasticity <- salt24h_both_DGE[(ratio <= 1.2 & ratio >=0.8),]
nrow(salt24h_S_D_plasticity)
#181 check the GO for those genes
salt24h_S_D_plasticity_high <- salt24h_S_D_plasticity[(abs(salt24h_S_D_plasticity$S_salt24h.log2FoldChange)>=4 & abs(salt24h_S_D_plasticity$D_salt24h.log2FoldChange)>=4),]
nrow(salt24h_S_D_plasticity_high)
#7

#make plot for the salt24h for 4 types of the genes
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt24h_DEG_S_D.pdf",width = 8,height = 8)
plot(salt24h_both_negative$S_salt24h.log2FoldChange,salt24h_both_negative$D_salt24h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-24h")
points(salt24h_both_positive$S_salt24h.log2FoldChange,salt24h_both_positive$D_salt24h.log2FoldChange,col="green")
points(salt24h_pS_nD$S_salt24h.log2FoldChange,salt24h_pS_nD$D_salt24h.log2FoldChange,col="blue")
points(salt24h_nS_pD$S_salt24h.log2FoldChange,salt24h_nS_pD$D_salt24h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#high
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt24h_DEG_S_D_high.pdf",width = 8,height = 8)
par(mar=c(5,4,4,2))
plot(salt24h_both_negative_high$S_salt24h.log2FoldChange,salt24h_both_negative_high$D_salt24h.log2FoldChange,col="orange",
     xlim=c(-25,25),ylim=c(-25,25),xlab="Sylcm",ylab="distcy",main="salt-24h")
points(salt24h_both_positive_high$S_salt24h.log2FoldChange,salt24h_both_positive_high$D_salt24h.log2FoldChange,col="green")
points(salt24h_pS_nD_high$S_salt24h.log2FoldChange,salt24h_pS_nD_high$D_salt24h.log2FoldChange,col="blue")
points(salt24h_nS_pD_high$S_salt24h.log2FoldChange,salt24h_nS_pD_high$D_salt24h.log2FoldChange,col="purple")
abline(h=0,lty=3,col="grey")
abline(v=0,lty=3,col="grey")
dev.off()
#3233


salt24h_S_unique <- salt24h[(salt24h$S_salt24h.padj<0.05 & salt24h$D_salt24h.padj>0.05),]
nrow(salt24h_S_unique)
#1472
salt24h_S_unique_high <- salt24h_S_unique[(abs(salt24h_S_unique$S_salt24h.log2FoldChange)>=4),]
nrow(salt24h_S_unique_high)
#17
salt24h_D_unique <- salt24h[(salt24h$S_salt24h.padj>0.05 & salt24h$D_salt24h.padj<0.05),]
nrow(salt24h_D_unique)
#4996
salt24h_D_unique_high <- salt24h_D_unique[(abs(salt24h_D_unique$D_salt24h.log2FoldChange>=4)),]
nrow(salt24h_D_unique_high)
#356

salt24h_both_NOT_DGE <- salt24h[(salt24h$S_salt24h.padj>0.05 & salt24h$D_salt24h.padj>0.05),]
nrow(salt24h_both_NOT_DGE)
#8833

#plot the barplot for the different categories for salt stress;

# Grouped Bar Plot
salt_plasticity <- c(nrow(salt1h_S_D_plasticity)/nrow(salt1h),
                     nrow(salt2h_S_D_plasticity)/nrow(salt2h),
                     nrow(salt5h_S_D_plasticity)/nrow(salt5h),
                     nrow(salt10h_S_D_plasticity)/nrow(salt10h),
                     nrow(salt24h_S_D_plasticity)/nrow(salt24h)
)

#high
salt_plasticity_high <- c(nrow(salt1h_S_D_plasticity_high)/nrow(salt1h),
                          nrow(salt2h_S_D_plasticity_high)/nrow(salt2h),
                          nrow(salt5h_S_D_plasticity_high)/nrow(salt5h),
                          nrow(salt10h_S_D_plasticity_high)/nrow(salt10h),
                          nrow(salt24h_S_D_plasticity_high)/nrow(salt24h)
)

salt_same_direction <- c(((nrow(salt1h_both_negative)+nrow(salt1h_both_positive)))/nrow(salt1h),
                         ((nrow(salt2h_both_negative)+nrow(salt2h_both_positive)))/nrow(salt2h),
                         ((nrow(salt5h_both_negative)+nrow(salt5h_both_positive)))/nrow(salt5h),
                         ((nrow(salt10h_both_negative)+nrow(salt10h_both_positive)))/nrow(salt10h),
                         ((nrow(salt24h_both_negative)+nrow(salt24h_both_positive)))/nrow(salt24h)
)

#high
salt_same_direction_high <- c(((nrow(salt1h_both_negative_high)+nrow(salt1h_both_positive_high)))/nrow(salt1h),
                              ((nrow(salt2h_both_negative_high)+nrow(salt2h_both_positive_high)))/nrow(salt2h),
                              ((nrow(salt5h_both_negative_high)+nrow(salt5h_both_positive_high)))/nrow(salt5h),
                              ((nrow(salt10h_both_negative_high)+nrow(salt10h_both_positive_high)))/nrow(salt10h),
                              ((nrow(salt24h_both_negative_high)+nrow(salt24h_both_positive_high)))/nrow(salt24h)
)

salt_opposite_direction <- c(((nrow(salt1h_pS_nD)+nrow(salt1h_nS_pD)))/nrow(salt1h),
                             ((nrow(salt2h_pS_nD)+nrow(salt2h_nS_pD)))/nrow(salt2h),
                             ((nrow(salt5h_pS_nD)+nrow(salt5h_nS_pD)))/nrow(salt5h),
                             ((nrow(salt10h_pS_nD)+nrow(salt10h_nS_pD)))/nrow(salt10h),
                             ((nrow(salt24h_pS_nD)+nrow(salt24h_nS_pD)))/nrow(salt24h)
)
#high
salt_opposite_direction_high <- c(((nrow(salt1h_pS_nD_high)+nrow(salt1h_nS_pD_high)))/nrow(salt1h),
                                  ((nrow(salt2h_pS_nD_high)+nrow(salt2h_nS_pD_high)))/nrow(salt2h),
                                  ((nrow(salt5h_pS_nD_high)+nrow(salt5h_nS_pD_high)))/nrow(salt5h),
                                  ((nrow(salt10h_pS_nD_high)+nrow(salt10h_nS_pD_high)))/nrow(salt10h),
                                  ((nrow(salt24h_pS_nD_high)+nrow(salt24h_nS_pD_high)))/nrow(salt24h)
)

salt_S_unique <- c(nrow(salt1h_S_unique)/nrow(salt1h),
                   nrow(salt2h_S_unique)/nrow(salt2h),
                   nrow(salt5h_S_unique)/nrow(salt5h),
                   nrow(salt10h_S_unique)/nrow(salt10h),
                   nrow(salt24h_S_unique)/nrow(salt24h)
)

#high
salt_S_unique_high <- c(nrow(salt1h_S_unique_high)/nrow(salt1h),
                        nrow(salt2h_S_unique_high)/nrow(salt2h),
                        nrow(salt5h_S_unique_high)/nrow(salt5h),
                        nrow(salt10h_S_unique_high)/nrow(salt10h),
                        nrow(salt24h_S_unique_high)/nrow(salt24h)
)

salt_D_unique <- c(nrow(salt1h_D_unique)/nrow(salt1h),
                   nrow(salt2h_S_unique)/nrow(salt2h),
                   nrow(salt5h_S_unique)/nrow(salt5h),
                   nrow(salt10h_S_unique)/nrow(salt10h),
                   nrow(salt24h_S_unique)/nrow(salt24h)
)
#high
salt_D_unique_high <- c(nrow(salt1h_D_unique_high)/nrow(salt1h),
                        nrow(salt2h_S_unique_high)/nrow(salt2h),
                        nrow(salt5h_S_unique_high)/nrow(salt5h),
                        nrow(salt10h_S_unique_high)/nrow(salt10h),
                        nrow(salt24h_S_unique_high)/nrow(salt24h)
)

salt24h_both_NOT_DGE <- c(nrow(salt1h_both_NOT_DGE)/nrow(salt1h),
                          nrow(salt2h_both_NOT_DGE)/nrow(salt2h),
                          nrow(salt5h_both_NOT_DGE)/nrow(salt5h),
                          nrow(salt10h_both_NOT_DGE)/nrow(salt10h),
                          nrow(salt24h_both_NOT_DGE)/nrow(salt24h)
)

salt_catNb <- cbind(salt_plasticity,
                    salt_same_direction,
                    salt_opposite_direction,
                    salt_S_unique,
                    salt_D_unique,
                    salt24h_both_NOT_DGE)
#high
salt_catNb_high <- cbind(salt_plasticity_high,
                         salt_same_direction_high,
                         salt_opposite_direction_high,
                         salt_S_unique_high,
                         salt_D_unique_high)
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt_DEG_S_D.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(salt_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
barplot(salt_catNb,beside=T,col=(colfunc(5)),main="salt",ylim=c(0,0.8),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni","Both_exp_no_sig"))
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/salt_DEG_S_D_high.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(salt_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
barplot(salt_catNb_high,beside=T,col=(colfunc(5)),main="salt",ylim=c(0,0.03),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni"))
dev.off()

#legend
pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/legend.pdf",width = 8,height = 8)
colfunc<-colorRampPalette(c("yellow","red"))
#barplot(salt_catNb,beside=T,col=(colfunc(5)),legend = c("1h","2h","5h","10h","24h"))
par(mfrow=c(1, 1), mar=c(1, 1, 1, 1))
barplot(salt_catNb,beside=T,col=(colfunc(5)),main="salt",ylim=c(0,0.8),names.arg=c("Robust", "Same_dir", "Opp_dir","S_uni","D_uni","Both_exp_no_sig"))

legend("topleft", 
       legend = c("1h","2h","5h","10h","24h"), 
       fill = (colfunc(5)), ncol = 2,
       cex = 0.75)

dev.off()

###Find the union o intersection for each category, and do GO analysis
#Heat robust:
#take the overlaps 
robust_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_S_D_plasticity$Syl_ID,heat2h_S_D_plasticity$Syl_ID),heat5h_S_D_plasticity$Syl_ID),heat10h_S_D_plasticity$Syl_ID),heat24h_S_D_plasticity$Syl_ID)
#no overlapped among all five timepoints
robust_heat_union <- union(union(union(union(heat1h_S_D_plasticity$Syl_ID,heat2h_S_D_plasticity$Syl_ID),heat5h_S_D_plasticity$Syl_ID),heat10h_S_D_plasticity$Syl_ID),heat24h_S_D_plasticity$Syl_ID)
robust_heat_high_union <- union(union(union(union(heat1h_S_D_plasticity_high$Syl_ID,heat2h_S_D_plasticity_high$Syl_ID),heat5h_S_D_plasticity_high$Syl_ID),heat10h_S_D_plasticity_high$Syl_ID),heat24h_S_D_plasticity_high$Syl_ID)

#Heat both_negative
#take the overlaps 
both_negative_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_both_negative$Syl_ID,heat2h_both_negative$Syl_ID),heat5h_both_negative$Syl_ID),heat10h_both_negative$Syl_ID),heat24h_both_negative$Syl_ID)
#no overlapped among all five timepoints
both_negative_heat_union <- union(union(union(union(heat1h_both_negative$Syl_ID,heat2h_both_negative$Syl_ID),heat5h_both_negative$Syl_ID),heat10h_both_negative$Syl_ID),heat24h_both_negative$Syl_ID)
both_negative_heat_high_union <- union(union(union(union(heat1h_both_negative_high$Syl_ID,heat2h_both_negative_high$Syl_ID),heat5h_both_negative_high$Syl_ID),heat10h_both_negative_high$Syl_ID),heat24h_both_negative_high$Syl_ID)


#Heat both_positive
both_positive_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_both_positive$Syl_ID,heat2h_both_positive$Syl_ID),heat5h_both_positive$Syl_ID),heat10h_both_positive$Syl_ID),heat24h_both_positive$Syl_ID)
#no overlapped among all five timepoints
both_positive_heat_union <- union(union(union(union(heat1h_both_positive$Syl_ID,heat2h_both_positive$Syl_ID),heat5h_both_positive$Syl_ID),heat10h_both_positive$Syl_ID),heat24h_both_positive$Syl_ID)
both_positive_heat_union_high <- union(union(union(union(heat1h_both_positive_high$Syl_ID,heat2h_both_positive_high$Syl_ID),heat5h_both_positive_high$Syl_ID),heat10h_both_positive_high$Syl_ID),heat24h_both_positive_high$Syl_ID)

#Heat opposite direction
#Heat salt1h_pS_nD
pS_nD_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_pS_nD$Syl_ID,heat2h_pS_nD$Syl_ID),heat5h_pS_nD$Syl_ID),heat10h_pS_nD$Syl_ID),heat24h_pS_nD$Syl_ID)
#overlapped among all five timepoints
pS_nD_heat_union <- union(union(union(union(heat1h_pS_nD$Syl_ID,heat2h_pS_nD$Syl_ID),heat5h_pS_nD$Syl_ID),heat10h_pS_nD$Syl_ID),heat24h_pS_nD$Syl_ID)
pS_nD_heat_union_high <- union(union(union(union(heat1h_pS_nD_high$Syl_ID,heat2h_pS_nD_high$Syl_ID),heat5h_pS_nD_high$Syl_ID),heat10h_pS_nD_high$Syl_ID),heat24h_pS_nD_high$Syl_ID)

#Heat opposite direction
#Heat salt1h_nS_pD
nS_pD_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_nS_pD$Syl_ID,heat2h_nS_pD$Syl_ID),heat5h_nS_pD$Syl_ID),heat10h_nS_pD$Syl_ID),heat24h_nS_pD$Syl_ID)
#overlapped among all five timepoints
nS_pD_heat_union <- union(union(union(union(heat1h_nS_pD$Syl_ID,heat2h_nS_pD$Syl_ID),heat5h_nS_pD$Syl_ID),heat10h_nS_pD$Syl_ID),heat24h_nS_pD$Syl_ID)
nS_pD_heat_union_high <- union(union(union(union(heat1h_nS_pD_high$Syl_ID,heat2h_nS_pD_high$Syl_ID),heat5h_nS_pD_high$Syl_ID),heat10h_nS_pD_high$Syl_ID),heat24h_nS_pD_high$Syl_ID)

#S_unique heat1h_S_unique
S_unique_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_S_unique$Syl_ID,heat2h_S_unique$Syl_ID),heat5h_S_unique$Syl_ID),heat10h_S_unique$Syl_ID),heat24h_S_unique$Syl_ID)
S_unique_heat_intersect_high <- intersect(intersect(intersect(intersect(heat1h_S_unique_high$Syl_ID,heat2h_S_unique_high$Syl_ID),heat5h_S_unique_high$Syl_ID),heat10h_S_unique_high$Syl_ID),heat24h_S_unique_high$Syl_ID)

S_unique_heat_union <- union(union(union(union(heat1h_S_unique$Syl_ID,heat2h_S_unique$Syl_ID),heat5h_S_unique$Syl_ID),heat10h_S_unique$Syl_ID),heat24h_S_unique$Syl_ID)
S_unique_heat_union_high <- union(union(union(union(heat1h_S_unique_high$Syl_ID,heat2h_S_unique_high$Syl_ID),heat5h_S_unique_high$Syl_ID),heat10h_S_unique_high$Syl_ID),heat24h_S_unique_high$Syl_ID)

#D_unique heat1h_D_unique
D_unique_heat_intersect <- intersect(intersect(intersect(intersect(heat1h_D_unique$Syl_ID,heat2h_D_unique$Syl_ID),heat5h_D_unique$Syl_ID),heat10h_D_unique$Syl_ID),heat24h_D_unique$Syl_ID)
D_unique_heat_union <- union(union(union(union(heat1h_D_unique$Syl_ID,heat2h_D_unique$Syl_ID),heat5h_D_unique$Syl_ID),heat10h_D_unique$Syl_ID),heat24h_D_unique$Syl_ID)
D_unique_heat_union_high <- union(union(union(union(heat1h_D_unique_high$Syl_ID,heat2h_D_unique_high$Syl_ID),heat5h_D_unique_high$Syl_ID),heat10h_D_unique_high$Syl_ID),heat24h_D_unique_high$Syl_ID)

#Write files:
write.table(
  robust_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/robust_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  robust_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/robust_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_negative_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_negative_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_negative_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_negative_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_positive_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_positive_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_positive_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_positive_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  pS_nD_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pS_nD_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  pS_nD_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pS_nD_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  nS_pD_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/nS_pD_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  nS_pD_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/nS_pD_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_heat_union_high,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_heat_union_high.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  D_unique_heat_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/D_unique_heat_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  D_unique_heat_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/D_unique_heat_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

#drought!!!!!!!!!!
#drought robust:
#take the overlaps 
robust_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_S_D_plasticity$Syl_ID,drought2h_S_D_plasticity$Syl_ID),drought5h_S_D_plasticity$Syl_ID),drought10h_S_D_plasticity$Syl_ID),drought24h_S_D_plasticity$Syl_ID)
#no overlapped among all five timepoints
robust_drought_union <- union(union(union(union(drought1h_S_D_plasticity$Syl_ID,drought2h_S_D_plasticity$Syl_ID),drought5h_S_D_plasticity$Syl_ID),drought10h_S_D_plasticity$Syl_ID),drought24h_S_D_plasticity$Syl_ID)
#drought both_negative
#take the overlaps 
both_negative_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_both_negative$Syl_ID,drought2h_both_negative$Syl_ID),drought5h_both_negative$Syl_ID),drought10h_both_negative$Syl_ID),drought24h_both_negative$Syl_ID)
#no overlapped among all five timepoints
both_negative_drought_union <- union(union(union(union(drought1h_both_negative$Syl_ID,drought2h_both_negative$Syl_ID),drought5h_both_negative$Syl_ID),drought10h_S_D_plasticity$Syl_ID),drought24h_S_D_plasticity$Syl_ID)
#drought both_positive
both_positive_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_both_positive$Syl_ID,drought2h_both_positive$Syl_ID),drought5h_both_positive$Syl_ID),drought10h_both_positive$Syl_ID),drought24h_both_positive$Syl_ID)
#no overlapped among all five timepoints
both_positive_drought_union <- union(union(union(union(drought1h_both_positive$Syl_ID,drought2h_both_positive$Syl_ID),drought5h_both_positive$Syl_ID),drought10h_both_positive$Syl_ID),drought24h_both_positive$Syl_ID)
#drought opposite direction
#drought salt1h_pS_nD
pS_nD_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_pS_nD$Syl_ID,drought2h_pS_nD$Syl_ID),drought5h_pS_nD$Syl_ID),drought10h_pS_nD$Syl_ID),drought24h_pS_nD$Syl_ID)
#overlapped among all five timepoints
pS_nD_drought_union <- union(union(union(union(drought1h_pS_nD$Syl_ID,drought2h_pS_nD$Syl_ID),drought5h_pS_nD$Syl_ID),drought10h_pS_nD$Syl_ID),drought24h_pS_nD$Syl_ID)

#drought opposite direction
#drought salt1h_nS_pD
nS_pD_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_nS_pD$Syl_ID,drought2h_nS_pD$Syl_ID),drought5h_nS_pD$Syl_ID),drought10h_nS_pD$Syl_ID),drought24h_nS_pD$Syl_ID)
#overlapped among all five timepoints
nS_pD_drought_union <- union(union(union(union(drought1h_nS_pD$Syl_ID,drought2h_nS_pD$Syl_ID),drought5h_nS_pD$Syl_ID),drought10h_nS_pD$Syl_ID),drought24h_nS_pD$Syl_ID)

#S_unique drought1h_S_unique
S_unique_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_S_unique$Syl_ID,drought2h_S_unique$Syl_ID),drought5h_S_unique$Syl_ID),drought10h_S_unique$Syl_ID),drought24h_S_unique$Syl_ID)
S_unique_drought_union <- union(union(union(union(drought1h_S_unique$Syl_ID,drought2h_S_unique$Syl_ID),drought5h_S_unique$Syl_ID),drought10h_S_unique$Syl_ID),drought24h_S_unique$Syl_ID)

#D_unique drought1h_D_unique
D_unique_drought_intersect <- intersect(intersect(intersect(intersect(drought1h_D_unique$Syl_ID,drought2h_D_unique$Syl_ID),drought5h_D_unique$Syl_ID),drought10h_D_unique$Syl_ID),drought24h_D_unique$Syl_ID)
D_unique_drought_union <- union(union(union(union(drought1h_D_unique$Syl_ID,drought2h_D_unique$Syl_ID),drought5h_D_unique$Syl_ID),drought10h_D_unique$Syl_ID),drought24h_D_unique$Syl_ID)

#Write files:
write.table(
  robust_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/robust_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  robust_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/robust_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_negative_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_negative_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_negative_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_negative_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_positive_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_positive_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_positive_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_positive_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  pS_nD_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pS_nD_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  pS_nD_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pS_nD_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  nS_pD_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/nS_pD_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  nS_pD_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/nS_pD_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  D_unique_drought_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/D_unique_drought_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  D_unique_drought_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/D_unique_drought_union.txt",
  quote=F,
  row.names=T,
  sep="\t")
#Salt
#salt robust:
#take the overlaps 
robust_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_S_D_plasticity$Syl_ID,salt2h_S_D_plasticity$Syl_ID),salt5h_S_D_plasticity$Syl_ID),salt10h_S_D_plasticity$Syl_ID),salt24h_S_D_plasticity$Syl_ID)
#no overlapped among all five timepoints
robust_salt_union <- union(union(union(union(salt1h_S_D_plasticity$Syl_ID,salt2h_S_D_plasticity$Syl_ID),salt5h_S_D_plasticity$Syl_ID),salt10h_S_D_plasticity$Syl_ID),salt24h_S_D_plasticity$Syl_ID)
#salt both_negative
#take the overlaps 
both_negative_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_both_negative$Syl_ID,salt2h_both_negative$Syl_ID),salt5h_both_negative$Syl_ID),salt10h_both_negative$Syl_ID),salt24h_both_negative$Syl_ID)
#no overlapped among all five timepoints
both_negative_salt_union <- union(union(union(union(salt1h_both_negative$Syl_ID,salt2h_both_negative$Syl_ID),salt5h_both_negative$Syl_ID),salt10h_S_D_plasticity$Syl_ID),salt24h_S_D_plasticity$Syl_ID)
#salt both_positive
both_positive_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_both_positive$Syl_ID,salt2h_both_positive$Syl_ID),salt5h_both_positive$Syl_ID),salt10h_both_positive$Syl_ID),salt24h_both_positive$Syl_ID)
#no overlapped among all five timepoints
both_positive_salt_union <- union(union(union(union(salt1h_both_positive$Syl_ID,salt2h_both_positive$Syl_ID),salt5h_both_positive$Syl_ID),salt10h_both_positive$Syl_ID),salt24h_both_positive$Syl_ID)
#salt opposite direction
#salt salt1h_pS_nD
pS_nD_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_pS_nD$Syl_ID,salt2h_pS_nD$Syl_ID),salt5h_pS_nD$Syl_ID),salt10h_pS_nD$Syl_ID),salt24h_pS_nD$Syl_ID)
#overlapped among all five timepoints
pS_nD_salt_union <- union(union(union(union(salt1h_pS_nD$Syl_ID,salt2h_pS_nD$Syl_ID),salt5h_pS_nD$Syl_ID),salt10h_pS_nD$Syl_ID),salt24h_pS_nD$Syl_ID)

#salt opposite direction
#salt salt1h_nS_pD
nS_pD_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_nS_pD$Syl_ID,salt2h_nS_pD$Syl_ID),salt5h_nS_pD$Syl_ID),salt10h_nS_pD$Syl_ID),salt24h_nS_pD$Syl_ID)
#overlapped among all five timepoints
nS_pD_salt_union <- union(union(union(union(salt1h_nS_pD$Syl_ID,salt2h_nS_pD$Syl_ID),salt5h_nS_pD$Syl_ID),salt10h_nS_pD$Syl_ID),salt24h_nS_pD$Syl_ID)

#S_unique salt1h_S_unique
S_unique_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_S_unique$Syl_ID,salt2h_S_unique$Syl_ID),salt5h_S_unique$Syl_ID),salt10h_S_unique$Syl_ID),salt24h_S_unique$Syl_ID)
S_unique_salt_union <- union(union(union(union(salt1h_S_unique$Syl_ID,salt2h_S_unique$Syl_ID),salt5h_S_unique$Syl_ID),salt10h_S_unique$Syl_ID),salt24h_S_unique$Syl_ID)

#D_unique salt1h_D_unique
D_unique_salt_intersect <- intersect(intersect(intersect(intersect(salt1h_D_unique$Syl_ID,salt2h_D_unique$Syl_ID),salt5h_D_unique$Syl_ID),salt10h_D_unique$Syl_ID),salt24h_D_unique$Syl_ID)
D_unique_salt_union <- union(union(union(union(salt1h_D_unique$Syl_ID,salt2h_D_unique$Syl_ID),salt5h_D_unique$Syl_ID),salt10h_D_unique$Syl_ID),salt24h_D_unique$Syl_ID)

#Write files:
write.table(
  robust_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/robust_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  robust_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/robust_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_negative_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_negative_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_negative_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_negative_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_positive_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_positive_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  both_positive_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/both_positive_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  pS_nD_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pS_nD_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  pS_nD_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/pS_nD_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  nS_pD_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/nS_pD_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  nS_pD_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/nS_pD_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  S_unique_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/S_unique_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  D_unique_salt_intersect,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/D_unique_salt_intersect.txt",
  quote=F,
  row.names=T,
  sep="\t")

write.table(
  D_unique_salt_union,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/D_unique_salt_union.txt",
  quote=F,
  row.names=T,
  sep="\t")

###

#quantile of the foldchanges for the treatment vs ctr
quantile(heat1h_both_DGE$S_heat1h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-21.8319699  -2.8078305  -2.1015597  -1.6470801  -1.3476695  -1.1029811  -0.9151323 
#35%         40%         45%         50%         55%         60%         65% 
#-0.7415979  -0.5884648  -0.4419605   0.3036084   0.4701260   0.6057945   0.7534609 
#70%         75%         80%         85%         90%         95%        100% 
#0.9279660   1.1502611   1.4487836   1.9106211   2.6647741   4.0031428  16.8373368 
quantile(heat1h_both_DGE$D_heat1h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-24.5451733  -2.8616761  -2.2077706  -1.8107325  -1.4632684  -1.2183770  -1.0186985 
#35%         40%         45%         50%         55%         60%         65% 
#-0.8279067  -0.6701950  -0.5002511  -0.2311667   0.4727587   0.6103410   0.7696414 
#70%         75%         80%         85%         90%         95%        100% 
#0.9497638   1.2116362   1.5907121   2.0601683   2.9340040   4.3307349  14.8090518
quantile(heat2h_both_DGE$S_heat2h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-24.0369596  -2.8553700  -2.0876465  -1.6660972  -1.3413175  -1.0934228  -0.9202391 
#35%         40%         45%         50%         55%         60%         65% 
#-0.7689992  -0.6190208  -0.4948132  -0.3367292   0.3869421   0.5114587   0.6179104 
#70%         75%         80%         85%         90%         95%        100% 
#0.7375774   0.8779930   1.0570939   1.2979675   1.7066243   2.6583201  19.8774742 
quantile(heat2h_both_DGE$D_heat2h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-23.4193960  -2.9019312  -2.1618164  -1.7532620  -1.4381391  -1.1805808  -0.9732742 
#35%         40%         45%         50%         55%         60%         65% 
#-0.8185141  -0.6288902  -0.4983990  -0.3081928   0.4238086   0.5471699   0.6609463 
#70%         75%         80%         85%         90%         95%        100% 
#0.7941007   0.9484173   1.1441400   1.4564553   1.9398127   2.9546610  15.9084092 
quantile(heat5h_both_DGE$S_heat5h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-16.1340315  -2.0753923  -1.5203825  -1.2176088  -1.0101864  -0.8713318  -0.7340745 
#35%         40%         45%         50%         55%         60%         65% 
#-0.6152805  -0.4987720  -0.3390654   0.3487648   0.4567943   0.5426487   0.6303189 
#70%         75%         80%         85%         90%         95%        100% 
#0.7245163   0.8509785   0.9796096   1.1990060   1.5405948   2.4121010  12.1310346 
quantile(heat5h_both_DGE$D_heat5h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-16.5559031  -2.2009049  -1.5895621  -1.2812686  -1.0581670  -0.8666312  -0.7513187 
#35%         40%         45%         50%         55%         60%         65% 
#-0.6281397  -0.5143315  -0.3694153   0.3709092   0.4873768   0.5632849   0.6521008 
#70%         75%         80%         85%         90%         95%        100% 
#0.7685035   0.9048272   1.0866301   1.3085464   1.6893772   2.5015781  11.8960133 
quantile(heat10h_both_DGE$S_heat10h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-11.0663294  -2.0188028  -1.5307438  -1.2362853  -1.0378600  -0.8734200  -0.7467582 
#35%         40%         45%         50%         55%         60%         65% 
#-0.6072982  -0.4745218  -0.3518763   0.3491381   0.4493445   0.5308949   0.6290082 
#70%         75%         80%         85%         90%         95%        100% 
#0.7259118   0.8582663   1.0696016   1.3303397   1.8124028   2.9811939  10.3174755
quantile(heat10h_both_DGE$D_heat10h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -7.3220211 -2.1767962 -1.5628119 -1.2298782 -1.0344091 -0.8584568 -0.7119281 -0.5795860 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  -0.4330878  0.2897682  0.4264754  0.5146964  0.6026847  0.7039138  0.8432514  1.0096459 
#80%        85%        90%        95%       100% 
#1.2190982  1.5330052  1.9726173  3.0090624 22.7365448 
quantile(heat24h_both_DGE$S_heat24h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-24.5666040  -2.2463470  -1.6387118  -1.3421000  -1.1314820  -0.9380035  -0.7850869 
#35%         40%         45%         50%         55%         60%         65% 
#-0.6398805  -0.5096640  -0.3914918   0.3423894   0.4959190   0.6247014   0.7583252 
#70%         75%         80%         85%         90%         95%        100% 
#0.9195521   1.0821697   1.2967707   1.6085797   2.1722314   3.0946988  22.8750650 
quantile(heat24h_both_DGE$D_heat24h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-27.7770731  -2.3005854  -1.6497841  -1.3549799  -1.1185797  -0.9277049  -0.7790199 
#35%         40%         45%         50%         55%         60%         65% 
#-0.6224366  -0.4735705  -0.3101696   0.4470427   0.5942699   0.7122803   0.8365854 
#70%         75%         80%         85%         90%         95%        100% 
#1.0016928   1.1818925   1.4170390   1.7048731   2.1899530   3.2126889  17.9727398 


###Drought
#quantile of the foldchanges for the treatment vs ctr
quantile(drought1h_both_DGE$S_drought1h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -5.2095928 -2.1061400 -1.4597727 -1.0886332 -0.8290386 -0.6184087 -0.4511838  0.3520069 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.4911045  0.6186544  0.7765164  0.9560282  1.2025599  1.4334623  1.7062001  2.0866233 
#80%        85%        90%        95%       100% 
#2.5311882  3.1222292  3.9702414  5.3522791 10.7592596
quantile(drought1h_both_DGE$D_drought1h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-10.6608684  -3.1713495  -2.4822718  -2.0031798  -1.6872418  -1.4775051  -1.2853541 
#35%         40%         45%         50%         55%         60%         65% 
#-1.1010483  -0.8749657   0.5924687   1.0144743   1.2969561   1.5232501   1.8141775 
#70%         75%         80%         85%         90%         95%        100% 
#2.0868485   2.3763543   2.8079635   3.2069278   3.9298836   5.3766373  23.1697311 
quantile(drought2h_both_DGE$S_drought2h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -6.6301585 -2.5721603 -1.7737382 -1.3229960 -0.9520445 -0.7580313 -0.5451498  0.2754221 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.6394406  0.8992142  1.1783918  1.4797786  1.7660204  2.0818065  2.4980446  2.9600853 
#80%        85%        90%        95%       100% 
#3.4685938  4.1754960  5.2255157  6.7758739 27.5773384 
quantile(drought2h_both_DGE$D_drought2h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-23.4193960  -2.9019312  -2.1618164  -1.7532620  -1.4381391  -1.1805808  -0.9732742 
#35%         40%         45%         50%         55%         60%         65% 
#-0.8185141  -0.6288902  -0.4983990  -0.3081928   0.4238086   0.5471699   0.6609463 
#70%         75%         80%         85%         90%         95%        100% 
#0.7941007   0.9484173   1.1441400   1.4564553   1.9398127   2.9546610  15.9084092 
quantile(drought5h_both_DGE$S_drought5h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -8.1835310 -2.8195953 -2.3586055 -1.9655275 -1.7047736 -1.4617084 -1.2358841 -1.0022062 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.7492114  1.2032362  1.4334867  1.6933688  1.9562312  2.2011837  2.4823091  2.8052389 
#80%        85%        90%        95%       100%  
quantile(drought5h_both_DGE$D_drought5h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -8.8903194 -3.2816457 -2.5048717 -2.1308760 -1.8355619 -1.5996152 -1.3874261 -1.2125818 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  -1.0251165 -0.8170783  0.5496353  0.8273189  1.0145102  1.2008377  1.4086503  1.6566483 
#80%        85%        90%        95%       100% 
#1.9547293  2.3595067  2.9129106  4.1745807 27.9476290 
quantile(drought10h_both_DGE$S_drought10h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -7.2170881 -2.9346944 -2.1512317 -1.7051680 -1.3838679 -1.1226416 -0.8770172 -0.6679389 
#40%        45%        50%        55%        60%        65%        70%        75% 
##  -0.4386610  0.4407770  0.6885693  0.9078594  1.1644275  1.4291048  1.7406622  2.1142032 
#80%        85%        90%        95%       100% 
#2.5630894  3.2343079  4.1582140  5.8019271 26.6508088
quantile(drought10h_both_DGE$D_drought10h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -9.5851306 -3.4127163 -2.6509480 -2.2125978 -1.8769170 -1.6130530 -1.4008602 -1.1872372 
#40%        45%        50%        55%        60%        65%        70%        75% 
##  -0.9544750 -0.6629951  0.7477591  0.9726624  1.1682310  1.4044857  1.6445792  1.9446957 
#80%        85%        90%        95%       100% 
#2.2904696  2.7357268  3.4248142  4.8333428 29.0721013 
quantile(drought24h_both_DGE$S_drought24h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -8.6959955 -3.1987639 -2.3776917 -1.8984159 -1.4906240 -1.2004768 -0.9018863 -0.6316791 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.2788917  0.6617312  0.9131286  1.1915969  1.4439660  1.7928367  2.1789162  2.6227203 
#80%        85%        90%        95%       100% 
#3.2156935  3.9548017  4.9751037  6.8946399 27.3708547 
quantile(drought24h_both_DGE$D_drought24h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -5.6928397 -2.8626598 -2.2259505 -1.8871618 -1.6159126 -1.3913154 -1.1785979 -0.9578997 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  -0.6743716  0.7737352  1.0005719  1.2077041  1.4062789  1.6313572  1.8778080  2.1550018 
#80%        85%        90%        95%       100% 
#2.5414809  3.0921546  3.9206660  5.3003547 29.9152286
#quantile of the foldchanges for the treatment vs ctr

###salt
quantile(salt1h_both_DGE$S_salt1h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -3.6822636 -0.8682252 -0.6210929 -0.3992563  0.3979172  0.5208022  0.6358174  0.7493225 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.8249561  0.9355137  1.0300184  1.1691598  1.3037890  1.4576308  1.6306163  1.8468837 
#80%        85%        90%        95%       100% 
#2.1143286  2.4492384  2.9777195  3.5599389  7.8341858 
quantile(salt1h_both_DGE$D_salt1h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-18.2101049  -4.6260816  -3.6357140  -3.0769845  -2.7128434  -2.3949819  -2.1605941 
#35%         40%         45%         50%         55%         60%         65% 
#-1.9979646  -1.7866479  -1.6516212  -1.5195882  -1.3808978  -1.2411768  -1.0275077 
#70%         75%         80%         85%         90%         95%        100% 
#-0.8470234   0.7379449   1.1464112   1.5172330   1.8916107   2.4286871   5.8261848 
quantile(salt2h_both_DGE$S_salt2h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-26.3041303  -1.2440849  -0.8717260  -0.7411064  -0.6106500  -0.5416621  -0.4968636 
#35%         40%         45%         50%         55%         60%         65% 
#-0.4178523   0.3609338   0.4621003   0.5907963   0.6910588   0.8352918   0.9902876 
#70%         75%         80%         85%         90%         95%        100% 
#1.1540063   1.3729306   1.6719433   2.0886296   2.5333799   3.3408418   5.7752664
quantile(salt2h_both_DGE$D_salt2h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -6.5814071 -3.1841936 -2.5450930 -2.1994309 -1.9366830 -1.7421553 -1.5640744 -1.4459384 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  -1.3113401 -1.1830502 -1.0980704 -1.0141892 -0.9202642 -0.8574281 -0.7430947 -0.5408977 
#80%        85%        90%        95%       100% 
#0.8290113  1.0929443  1.4632336  1.9727271 18.3218177 
quantile(salt5h_both_DGE$S_salt5h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%          5%         10%         15%         20%         25%         30% 
#-16.1396989  -0.9423253  -0.7503310  -0.6390733  -0.5406012  -0.4779672  -0.4473190 
#35%         40%         45%         50%         55%         60%         65% 
#-0.2953825   0.4244696   0.5199209   0.6278450   0.7213122   0.8076939   0.9417257 
#70%         75%         80%         85%         90%         95%        100% 
#1.0725045   1.1925606   1.3754539   1.4933845   1.8926889   2.3276729   4.7600163 
quantile(salt5h_both_DGE$D_salt5h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -4.4142174 -1.8807611 -1.5084052 -1.3115277 -1.2278180 -1.1049944 -1.0186921 -0.8641814 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  -0.7457679 -0.5574186  0.7757373  1.0549358  1.1551343  1.3868392  1.6091932  1.7607186 
#80%        85%        90%        95%       100% 
#1.9949700  2.2015916  2.5258456  3.2175554  7.8937835 
quantile(salt10h_both_DGE$S_salt10h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -5.5844026 -1.1548955 -0.7912254 -0.5721599 -0.4764596 -0.3612613  0.3738106  0.4367705 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.5116417  0.5501579  0.6241297  0.6995623  0.7818627  0.8771808  0.9602273  1.1513583 
#80%        85%        90%        95%       100% 
#1.3704737  1.6847574  2.3318398  3.2875181 18.9392082
quantile(salt10h_both_DGE$D_salt10h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -7.4007871 -2.3908327 -1.9429300 -1.6871030 -1.4729634 -1.2566188 -1.0250191 -0.8587919 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.5445933  0.7544387  0.8915841  1.0477949  1.1804072  1.3612561  1.5540496  1.7860474 
#80%        85%        90%        95%       100% 
#2.0593292  2.5254488  3.1674082  4.1473615 23.0366376 
quantile(salt24h_both_DGE$S_salt24h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -7.4007871 -2.3908327 -1.9429300 -1.6871030 -1.4729634 -1.2566188 -1.0250191 -0.8587919 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  0.5445933  0.7544387  0.8915841  1.0477949  1.1804072  1.3612561  1.5540496  1.7860474 
#80%        85%        90%        95%       100% 
#2.0593292  2.5254488  3.1674082  4.1473615 23.0366376  
quantile(salt24h_both_DGE$D_salt24h.log2FoldChange, probs = seq(0, 1, 0.05))
#0%         5%        10%        15%        20%        25%        30%        35% 
#  -6.3079028 -2.5865958 -2.2420317 -2.0022513 -1.8308983 -1.6618968 -1.5068045 -1.3626047 
#40%        45%        50%        55%        60%        65%        70%        75% 
#  -1.2294118 -1.0979235 -0.9281924 -0.7440035  0.6673355  1.0282003  1.2719411  1.6307124 
#80%        85%        90%        95%       100% 
#1.9285661  2.4613111  3.2749736  4.5643948 26.1641040 
##Based on the quantile of the foldchanges, I would just pull out the foldchabges >=4 genes and see each categories they fall down:

