setwd("D:/Rfiles/tetfig4")
library(vegan)
library(ggplot2)
library(reshape2)
library(gridExtra)

##barplot
#数据读取和整理
sdr <- read.csv("sdr.csv", row.names = 1, header = T, sep = ",")#DNA和RNA中检测的Sepceis
pdr <- read.csv("pdr.csv", row.names = 1, header = T, sep = ",")#DNA和RNA中检测的Pathway
g <- read.csv("gs.csv", header = T, sep = ",")
g$group <- factor(g$group, levels = c("H1D","H2D","H3D","H4D","H5D","H1R","H2R","H3R","H4R","H5R"))
dm <- function(data){
  data$mean <- rowMeans(data)
  temp <- data[order(data$mean, decreasing= T),]
  data_major <- temp[1:20,-ncol(data)]
  data_major <- data.frame(t(data_major))
  data_major$Others <- 100 - rowSums(data_major)
  data_major <- data.frame(t(data_major))
  return(data_major)
}##提取相对丰度前20的物种或通路
sdr <- dm(sdr)
pdr <- dm(pdr)
sdr$Species <- factor(rownames(sdr), levels = rev(rownames(sdr)))
pdr$Pathway <- factor(rownames(pdr), levels = rev(rownames(pdr)))

sdrl <- melt(sdr)
sdrl$idx <- 1:nrow(sdrl)
names(g)[1] <- "variable"
sdrl <- merge(sdrl, g, by = 'variable')
sdrl <- sdrl[order(sdrl$idx),]
sdrl <- sdrl[,-4]

pdrl <- melt(pdr)
pdrl$idx <- 1:nrow(pdrl)
names(g)[1] <- "variable"
pdrl <- merge(pdrl, g, by = 'variable')
pdrl <- pdrl[order(pdrl$idx),]
pdrl <- pdrl[,-4]
cols <- c("#447B66","#96CCA8","#C0EEDA","#CBD589","#86A845","#B3C648",
          "#FED46E","#FFDF6F","#F4D160","#8AC4D0","#58A3BC","#3E83A8",
          "#28527A","#134080","#1D6590","#A8D4E0","#537496","#55A6CB",
          "#A378B5","#766092","#DAB3DA","#D986B1","#F2A6C2","#FBC1AD",
          "#EF978F","#EFC2AD","#F56E4A","#F56E60","#F56E78","#D9D9D9")

#画图
p <- ggplot(sdrl, aes(variable, 100 * value, fill = Species)) +
  geom_col(position = 'stack', width = 0.6) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 12),legend.key.size=unit(0.2,'cm'))+
  scale_fill_manual(values =  rev(cols)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+facet_wrap(~group, scales = 'free_x', ncol = 5) +
  theme(strip.text = element_text(size = 12))

p1 <- ggplot(pdrl, aes(variable, 100 * value, fill = Pathway)) +
  geom_col(position = 'stack', width = 0.6) +
  labs(x = '', y = 'Relative Abundance(%)') +
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 13)) +
  theme(legend.text = element_text(size = 12, ncol(1)),legend.key.size=unit(0.2,'cm'))+
  scale_fill_manual(values =  rev(cols)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())+facet_wrap(~group, scales = 'free_x', ncol = 5) +
  theme(strip.text = element_text(size = 12))
# ggsave('pathway.pdf', p1, width = 16, height = 10, dpi = 1080)

##Bray-Curtis to FB
plotstackzx <- function(data1,data2,data3,data4,data5){
  cols <- c("#7DABD0","#CFE7EA","#FBC1AD","#F56E4A","#FFDF6F")
  time <- c("PreF","PostI24h","PostI36h","PostI48h","PostI60h","PostI72h","PostI84h",
            "PostI96h","PostI108h","PostI120h","PostI132h","PostI144h","PostI156h")
  t <- c(0,24,36,48,60,72,84,96,108,120,132,144,156)
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
  index2 <- t[1:len2]
  dist2 <- c(0,d2[1:(len2-1)])
  group2 <- rep('H2',each = len2)
  bcd2 <- data.frame(id=id2, group=group2, index=index2, dist=dist2)
  
  d3 <- vegdist(data3[,-1],"bray")
  len3 <- nrow(data3)
  id3 <- time[1:len3]
  index3 <- t[1:len3]
  dist3 <- c(0,d3[1:(len3-1)])
  group3 <- rep('H3',each = len3)
  bcd3 <- data.frame(id=id3, group=group3, index=index3, dist=dist3)
  
  d4 <- vegdist(data4[,-1],"bray")
  len4 <- nrow(data4)
  id4 <- time[1:len4]
  index4 <- t[1:len4]
  dist4 <- c(0,d4[1:(len4-1)])
  group4 <- rep('H4',each = len4)
  bcd4 <- data.frame(id=id4, group=group4, index=index4, dist=dist4)
  
  d5 <- vegdist(data5[,-1],"bray")
  len5 <- nrow(data5)
  id5 <- time[1:len5]
  index5 <- t[1:len5]
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
}#绘制折线图

#Species--ilealcecum--metagenome
simg <- read.csv("D:/Rfiles/simg.csv",header=T,sep=",")
dim(simg)
la <- nrow(simg)+1
n <- c(1)
j = 1
for (i in 1:nrow(simg)){
  if(simg[i,1] ==0){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
hc007ft <- data.frame(simg[n[1]:(n[2]-1),])
hc008ft <- data.frame(simg[n[2]:(n[3]-1),])
hc015ft <- data.frame(simg[n[3]:(n[4]-1),])
hc016ft <- data.frame(simg[n[4]:(n[5]-1),])
hc017ft <- data.frame(simg[n[5]:(n[6]-1),])
p1 <- plotstackzx(hc007simg,hc008simg,hc015simg,hc016simg,hc017simg)

#Species--ilealcecum--meta-transcriptome
simt <- read.csv("D:/Rfiles/simt.csv",header=T,sep=",")
dim(simt)
la <- nrow(simt)+1
n <- c(1)
j = 1
for (i in 1:nrow(simt)){
  if(simt[i,1] ==0){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
hc007simt <- data.frame(simt[n[1]:(n[2]-1),])
hc008simt <- data.frame(simt[n[2]:(n[3]-1),])
hc015simt <- data.frame(simt[n[3]:(n[4]-1),])
hc016simt <- data.frame(simt[n[4]:(n[5]-1),])
hc017simt <- data.frame(simt[n[5]:(n[6]-1),])
p2 <- plotstackzx(hc007simt,hc008simt,hc015simt,hc016simt,hc017simt)

#Pathway--ilealcecum--metagenome
pimg <- read.csv("D:/Rfiles/pimg.csv",header=T,sep=",")
dim(pimg)
la <- nrow(pimg)+1
n <- c(1)
j = 1
for (i in 1:nrow(pimg)){
  if(pimg[i,1] ==0){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
hc007pimg <- data.frame(pimg[n[1]:(n[2]-1),])
hc008pimg <- data.frame(pimg[n[2]:(n[3]-1),])
hc015pimg <- data.frame(pimg[n[3]:(n[4]-1),])
hc016pimg <- data.frame(pimg[n[4]:(n[5]-1),])
hc017pimg <- data.frame(pimg[n[5]:(n[6]-1),])
p3 <- plotstackzx(hc007pimg,hc008pimg,hc015pimg,hc016pimg,hc017pimg)

#Pathway--ilealcecum--meta-transcriptome
pimt <- read.csv("D:/Rfiles/pimt",header=T,sep=",")
dim(pimt)
la <- nrow(pimt)+1
n <- c(1)
j = 1
for (i in 1:nrow(pimt)){
  if(pimt[i,1] ==0){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
hc007pimt <- data.frame(pimt[n[1]:(n[2]-1),])
Hc008pimt <- data.frame(pimt[n[2]:(n[3]-1),])
hc015pimt <- data.frame(pimt[n[3]:(n[4]-1),])
hc016pimt <- data.frame(pimt[n[4]:(n[5]-1),])
hc017pimt <- data.frame(pimt[n[5]:(n[6]-1),])
p4 <- plotstackzx(hc007pimt,hc008pimt,hc015pimt,hc016pimt,hc017pimt)

pall <- grid.arrange(p1,p2,p3,p4,ncol = 1, nrow = 4) 
ggsave(paste("bray-curtis",".pdf", sep=""),pall,width = 6, height = 20,dpi = 1080)
