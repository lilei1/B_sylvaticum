 
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("topGO")
sessionInfo()
if (!requireNamespace("BiocManager", quietly=TRUE))
  +     install.packages("BiocManager")
BiocManager::install("Rgraphviz")


library("topGO")

geneID2GO <- readMappings(file = "/global/u2/l/llei2019/Brachypodium/Sylvaticum/TopGo/Bsylvaticum_490_v1.1_GO_gid.txt") 
geneUniverse <- names(geneID2GO) 
#geneUniverse
genesOfInterest <- read.delim("/global/u2/l/llei2019/Brachypodium/Sylvaticum/TopGo/target_GO",sep="\t",header=FALSE)
genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata 
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)
resultFisher <- runTest(myGOdata, algorithm="classic", statistic="fisher")
allRes <- GenTable(myGOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)


showSigOfNodes(myGOdata, score(resultFisher), firstSigNodes = 8, useInfo ='all')
output_file2 = paste("/global/u2/l/llei2019/Brachypodium/Sylvaticum/TopGo/","allDEG_TopGo", sep="_")
printGraph(myGOdata, resultFisher, firstSigNodes =8, fn.prefix = output_file2, useInfo = "all", pdfSW = TRUE)
length(usedGO(myGOdata))

  
myterms = c("GO:0019637", "GO:1901566", " GO:0051641","GO:0009987")
mygenes <- genesInTerm(myGOdata, myterms)
for (i in 1:length(myterms))
{
  myterm <- myterms[i]
  mygenesforterm <- mygenes[myterm][[1]]
  mygenesforterm <- paste(mygenesforterm, collapse=',')
  print(paste("Term",myterm,"genes:",mygenesforterm))
}
