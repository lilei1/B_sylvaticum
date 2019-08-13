install.packages("UpSetR")

library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")
head(movies)
nrow(movies)
upset(movies, sets = c("Action", "Adventure", "Comedy", "Drama", "Mystery", "Thriller", "Romance", "War"))


data <- read.delim(file="/Users/LiLei/Downloads/orthodbs_round2/sorted_masked_gene_ortholog.matrix", header = T, sep = "\t")
head(data)
pdf(file = "/Users/LiLei/Downloads/orthodbs_round2/upset_orth.pdf",width = 42,height = 15)
#upset(data,sets = c("Sbicolor","OsativaKitaake","PhalliiHAL","Bsylvaticum","Bstacei","Bdistachyon"), order.by = c("degree"), mainbar.y.label="Orthologous clusters in each intersection",sets.x.label
#="Orthologous clusters in each species",text.scale=2.5,point.size=7)

setsBarColors <- c('#4285F4','#4285F4', '#EA4335', '#FBBC05', '#FBBC05', '#FBBC05')

upset(data,
      nsets=6,
      nintersects=1000,
      sets = c("PhalliiHAL","Bsylvaticum","Sbicolor","Bstacei","Bdistachyon","OsativaKitaake"), 
      keep.order = TRUE, 
      mainbar.y.label="Orthologous clusters in each intersection",
      sets.x.label="Orthologous clusters in each species",
      text.scale=2.5,
      point.size=7,
      sets.bar.color=setsBarColors,
      queries = list(list(query = intersects,params = list("Bstacei","Bdistachyon","OsativaKitaake"), color ='#FBBC05' , active = T),
                     list(query = intersects,params = list("PhalliiHAL","Bsylvaticum"), color ='#4285F4' , active = T))
)

dev.off()

