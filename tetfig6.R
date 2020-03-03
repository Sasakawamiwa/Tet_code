setwd("D:/Rfiles/tetfig6")
mmajor <- read.csv("mmajor.csv", header = T, row.names = 1, sep = ",")
g <- read.csv("g.csv", header = T, sep = ",")
library(vegan)
library(ggplot2)
library(reshape2)
library(gridExtra)

##barplot
g$group <- factor(g$group, levels = c("H1F","H2F","H3F","H4F","H5F","H1U","H2U","H3U","H4U","H5U"))
mmajor$Metabolite <- factor(rownames(mmajor), levels = rev(rownames(mmajor)))

mmajor <- melt(mmajor)
mmajor$idx <- 1:nrow(mmajor)
names(g)[1] <- "variable"
mmajor <- merge(mmajor, g, by = 'variable')
mmajor <- mmajor[order(mmajor$idx),]
mmajor <- mmajor[,-4]

cols <- c("#447B66","#96CCA8","#C0EEDA","#CBD589","#86A845",
          "#FFDF6F","#F4D160","#F4D160","#58A3BC","#3E83A8",
          "#008B8B","#134080","#20B2AA","#A8D4E0","#537496","#55A6CB",
          "#A378B5","#766092","#DAB3DA","#D986B1","#F2A6C2","#FBC1AD",
          "#EF978F","#EFC2AD","#F56E4A","#F56E78","#D9D9D9")

p <- ggplot(mmajor, aes(variable, 100 * value, fill = Metabolite)) +
  geom_col(position = 'stack', width = 0.6) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 12),legend.key.size=unit(0.2,'cm'))+
  scale_fill_manual(values =  rev(cols)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+facet_wrap(~group, scales = 'free_x', ncol = 5) +
  theme(strip.text = element_text(size = 12))


##b-c
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
divert <- function(data,k){
  la <- nrow(data)+1
  n <- c(1)
  j = 1
  for (i in 1:nrow(data)){
    if(data[i,1] ==0){
      n[j] <- i
      j = j + 1
    }
  }
  n[j] <- la
  d <- data.frame(data[n[k]:(n[k+1]-1),])
}
plotstackzx <- function(data1,data2,data3,data4,data5){
  cols <- c("#7DABD0","#CFE7EA","#FBC1AD","#F56E4A","#FFDF6F")
  time <- c("PreF","PostF24h","PostF48h","PostF72h","PostF96h","PostF120h","PostF144h")
  d1 <- vegdist(data1[,-1],"bray")
  len1 <- nrow(data1)
  id1 <- time[1:len1]
  index1 <- data1[,1]
  dist1 <- c(0,d1[1:(len1-1)])
  group1 <- rep('H1',each = len1)
  bcd1 <- data.frame(id=id1, group=group1, index=index1, dist=dist1)
  
  d2 <- vegdist(data2[,-1],"bray")
  len2 <- nrow(data2)
  id2 <- time[1:len2]
  index2 <- data2[,1]
  dist2 <- c(0,d2[1:(len2-1)])
  group2 <- rep('H2',each = len2)
  bcd2 <- data.frame(id=id2, group=group2, index=index2, dist=dist2)
  
  d3 <- vegdist(data3[,-1],"bray")
  len3 <- nrow(data3)
  id3 <- time[1:len3]
  index3 <- data3[,1]
  dist3 <- c(0,d3[1:(len3-1)])
  group3 <- rep('H3',each = len3)
  bcd3 <- data.frame(id=id3, group=group3, index=index3, dist=dist3)
  
  d4 <- vegdist(data4[,-1],"bray")
  len4 <- nrow(data4)
  id4 <- time[1:len4]
  index4 <- data4[,1]
  dist4 <- c(0,d4[1:(len4-1)])
  group4 <- rep('H4',each = len4)
  bcd4 <- data.frame(id=id4, group=group4, index=index4, dist=dist4)
  
  d5 <- vegdist(data5[,-1],"bray")
  len5 <- nrow(data5)
  id5 <- time[1:len5]
  index5 <- data5[,1]
  dist5 <- c(0,d5[1:(len5-1)])
  group5 <- rep('H5',each = len5)
  bcd5 <- data.frame(id=id5, group=group5, index=index5, dist=dist5)
  
  db <- rbind(bcd1,bcd2,bcd3,bcd4,bcd5)
  ggplot(db, aes(x=index, y=dist, colour = group)) + geom_line() + geom_point(size=4, shape=20)+
    xlab("Sampling time") + ylab("The Bray–Curtis Dissimilarity index to the PreF sample")+
    theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),axis.text=element_text(angle = 90))+
    scale_x_continuous(breaks=db$index, labels = db$id)+
    scale_color_manual(values = cols)+
    theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'))
}

fm <- read.csv("fm.csv", header = T, sep = ",")
fm[,-1] <- normalization(fm[,-1])
fm1 <- divert(fm,1)
fm2 <- divert(fm,2)
fm3 <- divert(fm,3)
fm4 <- divert(fm,4)
fm5 <- divert(fm,5)
p1 <- plotstackzx(fm1,fm2,fm3,fm4,fm5)

um <- read.csv("um.csv", header = T, sep = ",")
um[,-1] <- normalization(um[,-1])
um1 <- divert(um,1)
um2 <- divert(um,2)
um3 <- divert(um,3)
um4 <- divert(um,4)
um5 <- divert(um,5)
p2 <- plotstackzx(um1,um2,um3,um4,um5)

pall <- grid.arrange(p1,p2,ncol = 1, nrow = 2) 
ggsave(paste("bray-curtis",".pdf", sep=""),pall,width = 6, height = 10,dpi = 1080)
