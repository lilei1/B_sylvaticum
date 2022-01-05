#Write by Li:
file <- read.delim(file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/sniffles.txt",sep="\t")
head(file)
#   Drop the unmapped chromosome
file <- file[(file$Chrom !="chrUn"),]
#Take a different categories:
#ins <- file[file$type=="INS",]
head(file)
#del <- file[file$type=="DEL",]
#dup <- file[file$type=="DUP",]
#inv <- file[file$type=="INV",]
#tra <- file[file$type=="TRA",]

file$len <- abs(file$stop2 - file$start)
new_file <- file[,c(1,2,6,11,21)]
head(new_file)
names(new_file)[1]<-"Chromosome"
names(new_file)[2]<-"LeftBP"
names(new_file)[3]<-"RightBP"
names(new_file)[4]<-"Type"
names(new_file)[5]<-"Length"
head(cmmb)
cmmb <-cmmb[order(cmmb$Chromosome, cmmb$LeftBP),]
ins <- cmmb[cmmb$Type=="INS",]
del <- cmmb[cmmb$Type=="DEL",]
dup <- cmmb[cmmb$Type=="DUP",]
inv <- cmmb[cmmb$Type=="INV",]
tra <- cmmb[cmmb$Type=="TRA",]

ins
del
dup
inv
tra
#   Save the new file!
write.table(
  cmmb,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/cmmb.txt",
  quote=F,
  row.names=F,
  sep="\t")

write.table(
  ins,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/ins.txt",
  quote=F,
  row.names=F,
  sep="\t")

write.table(
  del,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/del.txt",
  quote=F,
  row.names=F,
  sep="\t")

write.table(
  dup,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/dup.txt",
  quote=F,
  row.names=F,
  sep="\t")

write.table(
  inv,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/inv.txt",
  quote=F,
  row.names=F,
  sep="\t")

write.table(
  tra,
  file="/global/u2/l/llei2019/Brachypodium/Sylvaticum/RNAseq/tra.txt",
  quote=F,
  row.names=F,
  sep="\t")

