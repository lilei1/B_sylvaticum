library(ggplot2)
TTM <- read.delim("Dropbox (Personal)/Brachpodium/sylvaticum/TE/MS/Proportion/proportion_input_ggplot.txt",sep="\t",header=T)
head(TTM)
TTM$Species <- factor(TTM$Species,levels = TTM$Species)
library(RColorBrewer)
display.brewer.all()
#ggplot(data = TTM, aes(x = Species, y = Proportion, fill = Categories)) + geom_bar(stat="identity") + coord_flip()+scale_x_discrete(limits=TTM$Species)
stck <- ggplot(data = TTM, aes(x = factor(Species,levels=unique(Species)), y = Proportion, 
                       fill =factor(Categories, levels=unique(Categories)) )) + geom_bar(position = "fill",stat="identity", width=0.4) + coord_flip()+ xlab("") + ylab("") +theme_bw() +
  theme(axis.line.y= element_blank(),
        axis.line.x = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.text=element_text(size=rel(1.0)),
        axis.text.x = element_text(size = 20,face="bold"),axis.text.y = element_text(size = 20,face="bold.italic")) + guides(fill=guide_legend(title="")) 

  stck+scale_fill_manual(values=c("#00bfff", "#ffbf00","#ff8000","#ffff00","#bfff00","#00b3b3","#d966ff","#66ff33","#9933ff","#adad85","#0080ff","#8000ff","#ff00bf","#ff0040","#009900","#808000","#cc9900","#9999ff"))
  ggsave("/Users/LiLei/Dropbox (Personal)/Brachpodium/sylvaticum/TE/MS/Proportion/proportion.pdf", width = 20, height =12)
  
  ###
  
 MM <- read.delim("Dropbox (Personal)/Brachpodium/sylvaticum/TE/MS/Proportion/len_Proportion.txt",sep="\t",header=T)
  head(MM)
  #MM$Species <- factor(MM$Species,levels = MM$Species)
  #library(RColorBrewer)
  #display.brewer.all()
  ggplot(data = MM, aes(x = factor(Species,levels=unique(Species)), y = Length, fill =factor(Categories, levels=unique(Categories)))) + geom_bar(stat="identity") + coord_flip()
  
stc <- ggplot(data = MM, aes(x = factor(Species,levels=unique(Species)), y = Length, 
                                 fill =factor(Categories, levels=unique(Categories)) )) + 
  geom_bar(stat="identity", width=0.4) + coord_flip()+ xlab("") + ylab("") +theme_bw() +
    theme(axis.line.y= element_blank(),
          axis.line.x = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.text=element_text(size=rel(1.0)),
          axis.text.x = element_text(size = 20,face="bold"),axis.text.y = element_text(size = 20,face="bold.italic")) + guides(fill=guide_legend(title="")) 
  
  stc+scale_fill_manual(values=c("#00bfff", "#ffbf00","#ff8000","#ffff00","#bfff00","#00b3b3","#d966ff","#66ff33","#9933ff","#adad85","#0080ff","#8000ff","#ff00bf","#ff0040","#009900","#808000","#cc9900","#9999ff"))
  
  ggsave("/Users/LiLei/Dropbox (Personal)/Brachpodium/sylvaticum/TE/MS/Proportion/length.pdf", width = 20, height =12)
  