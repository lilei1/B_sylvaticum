
remotes::install_github("r-lib/rlang", build_vignettes = TRUE)
install.packages("UpSetR")
library(UpSetR)

data <- read.delim(file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/Abiostress_timepoint_DGE_sum_TF.txt", header = T, sep = "\t")
#data <- read.delim(file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/test.txt", header = T, sep = "\t")

head(data)

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_orth_all.pdf",width = 150,height = 15)
#upset(data,sets = c("cond_1v6_Sig", "cond_1v11_Sig", "cond_1v16_Sig", "cond_2v7_Sig", "cond_2v12_Sig", "cond_2v17_Sig", "cond_3v8_Sig", "cond_3v13_Sig", "cond_3v18_Sig", "cond_4v9_Sig", "cond_4v14_Sig","cond_4v19_Sig","cond_5v10_Sig", "cond_5v15_Sig", "cond_5v20_Sig"), order.by = c("degree"), mainbar.y.label="No. of DFG",sets.x.label="No. of DFG",text.scale=2.5,point.size=7)
#data <- as.data.frame(data)
upset(data,      
      nsets=6,
      nintersects=1000,
      sets = c("Salt_1h","Salt_2h","Salt_5h","Salt_10h","Salt_24h","Drought_1h","Drought_2h","Drought_5h","Drought_10h","Drought_24h","Heat_1h", "Heat_2h", "Heat_5h","Heat_10h","Heat_24h"), 
      keep.order = TRUE, 
      mainbar.y.label ="No. of DFG",
      sets.x.label ="No. of DFG in each treatment/timepoint",
      text.scale=1,
      point.size=2)
dev.off()


upset(data,sets = c("Salt_1h","Salt_2h","Salt_5h","Salt_10h","Salt_24h","Drought_1h","Drought_2h","Drought_5h","Drought_10h","Drought_24h","Heat_1h", "Heat_2h", "Heat_5h","Heat_10h","Heat_24h"), order.by = "freq",mainbar.y.label="No. of DFG",sets.x.label ="No. of DFG in each treatment/timepoint",text.scale=1,point.size=4)

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_Salt_time.pdf",width = 25,height = 15)
upset(data,sets = c("Salt_1h","Salt_2h","Salt_5h","Salt_10h","Salt_24h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each timepoint",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_Drought_time.pdf",width = 25,height = 15)
upset(data,sets = c("Drought_1h","Drought_2h","Drought_5h","Drought_10h","Drought_24h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each timepoint",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_heat_time.pdf",width = 25,height = 15)
upset(data,sets = c("Heat_1h", "Heat_2h", "Heat_5h","Heat_10h","Heat_24h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each timepoint",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_abiotic_1h.pdf",width = 25,height = 15)
upset(data,sets = c("Salt_1h", "Drought_1h", "Heat_1h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each treatment (1h)",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_abiotic_2h.pdf",width = 25,height = 15)
upset(data,sets = c("Salt_2h", "Drought_2h", "Heat_2h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each treatment (2h)",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_abiotic_5h.pdf",width = 25,height = 15)
upset(data,sets = c("Salt_5h", "Drought_5h", "Heat_5h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each treatment (5h)",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_abiotic_10h.pdf",width = 25,height = 15)
upset(data,sets = c("Salt_10h", "Drought_10h", "Heat_10h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each treatment (10h)",text.scale=2.5,point.size=6)
dev.off()

pdf(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/upset_abiotic_24h.pdf",width = 25,height = 15)
upset(data,sets = c("Salt_24h", "Drought_24h", "Heat_24h"), keep.order = TRUE, mainbar.y.label="No. of DGE",sets.x.label ="No. of DGE in each treatment (24h)",text.scale=2.5,point.size=6)
dev.off()



