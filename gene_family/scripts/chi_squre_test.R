mem_count <- read.delim("Dropbox (Personal)/Brachpodium/sylvaticum/gene_family/member_count.txt",sep="\t",header=T)
head(mem_count)
#Consider Sorghum as the perennial
mem_count$annual <- mem_count$Osativa + mem_count$Bdistachyon + mem_count$Bstacei
mem_count$perenial <- mem_count$Bsylvaticum + mem_count$Phallii + mem_count$Sbicolor
mem_count$annual_exp <- (mem_count$Osativa + mem_count$Bdistachyon + mem_count$Bstacei + mem_count$Bsylvaticum + mem_count$Phallii + mem_count$Sbicolor)/2
mem_count$perenial_exp <- (mem_count$Osativa + mem_count$Bdistachyon + mem_count$Bstacei + mem_count$Bsylvaticum + mem_count$Phallii + mem_count$Sbicolor)/2

test <- mem_count[1-30,]
head(test)


#results <- t(apply(test,1,function(x)
#  with(chisq.test(x[8:9],p=x[10:11]),c(statistic,p.value=p.value))))

mem_count$chisq     <- with(mem_count,(annual-annual_exp)^2/annual_exp + (perenial-perenial_exp)^2/perenial_exp)

mem_count$p.value   <- pchisq(mem_count$chisq,df=1, lower.tail=F)

write.table(mem_count, "Dropbox (Personal)/Brachpodium/sylvaticum/gene_family/member_count_chi2.txt",quote =F,col.names =T,row.names =F,sep="\t")

#Consider Sorghum as the annual:
mem_count <- read.delim("Dropbox (Personal)/Brachpodium/sylvaticum/gene_family/member_count.txt",sep="\t",header=T)
head(mem_count)
mem_count$annual <- mem_count$Osativa + mem_count$Bdistachyon + mem_count$Bstacei+ mem_count$Sbicolor
mem_count$perenial <- mem_count$Bsylvaticum + mem_count$Phallii 
mem_count$annual_exp <- (mem_count$Osativa + mem_count$Bdistachyon + mem_count$Bstacei + mem_count$Bsylvaticum + mem_count$Phallii + mem_count$Sbicolor)*4/6
mem_count$perenial_exp <- (mem_count$Osativa + mem_count$Bdistachyon + mem_count$Bstacei + mem_count$Bsylvaticum + mem_count$Phallii + mem_count$Sbicolor)*2/6
mem_count$chisq     <- with(mem_count,(annual-annual_exp)^2/annual_exp + (perenial-perenial_exp)^2/perenial_exp)

mem_count$p.value   <- pchisq(mem_count$chisq,df=1, lower.tail=F)

write.table(mem_count, "Dropbox (Personal)/Brachpodium/sylvaticum/gene_family/member_count_chi2_sorgum_an.txt",quote =F,col.names =T,row.names =F,sep="\t")

#only include Brachypodium
mem_count <- read.delim("Dropbox (Personal)/Brachpodium/sylvaticum/gene_family/member_count.txt",sep="\t",header=T)
head(mem_count)
mem_count$annual <- mem_count$Bdistachyon + mem_count$Bstacei
mem_count$perenial <- mem_count$Bsylvaticum
mem_count$annual_exp <- (mem_count$Bdistachyon + mem_count$Bstacei + mem_count$Bsylvaticum)*2/3
mem_count$perenial_exp <- (mem_count$Bdistachyon + mem_count$Bstacei + mem_count$Bsylvaticum)*1/3
mem_count$chisq     <- with(mem_count,(annual-annual_exp)^2/annual_exp + (perenial-perenial_exp)^2/perenial_exp)

mem_count$p.value   <- pchisq(mem_count$chisq,df=1, lower.tail=F)

write.table(mem_count, "Dropbox (Personal)/Brachpodium/sylvaticum/gene_family/member_count_chi2_brachy.txt",quote =F,col.names =T,row.names =F,sep="\t")
