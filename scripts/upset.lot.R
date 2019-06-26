install.packages("UpSetR")

library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")
head(movies)
nrow(movies)
upset(movies, sets = c("Action", "Adventure", "Comedy", "Drama", "Mystery", "Thriller", "Romance", "War"))


data <- read.delim(file="/Users/LLei/Downloads/orthodbs/sorted_masked_gene_ortholog.matrix", header = T, sep = "\t")
head(data)
pdf(file = "/Users/LLei/Downloads/orthodbs/upset_orth.pdf",width = 25,height = 15)
upset(data,sets = c("Bdistachyon", "Bstacei", "Bsylvaticum", "OsativaKitaake", "Sbicolor"), order.by = c("degree"), mainbar.y.label="Orthologous clusters in each intersection",sets.x.label
="Orthologous clusters in each species",text.scale=2.5,point.size=7)

dev.off()
