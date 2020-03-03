setwd("D:/Rfiles/tetfig5")

###Pathway
##计算阶段
pmg <- read.csv("D:/Rfiles/rhythmsg/pmg.csv", header = T, sep = ",")
la <- nrow(pmg)+1
n <- c(1)
j = 1
for (i in 1:nrow(pmg)){
  if(pmg[i,1] == 24){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
id1 <- data.frame(pmg[n[1]:(n[2]-1),])
id2 <- data.frame(pmg[n[2]:(n[3]-1),])
id3 <- data.frame(pmg[n[3]:(n[4]-1),])
id4 <- data.frame(pmg[n[4]:(n[5]-1),])
id5 <- data.frame(pmg[n[5]:(n[6]-1),])
sumzero <- function(x) sum(x == 0)
datafilter <- function(data){
  len <- nrow(data)
  l <- len/2 + 1
  data <- data[,c(1,which(apply(data[,-1],2,sumzero) < l)+1)]
  data <- t(data) 
  data <- data.frame(data)
  seq <- rep(c("+", "-"), times = 10)
  derta <- data[,-1] #derta <- id2[,-1]
  for (i in 1:(len-1)){
    derta[,i] <- data[,i+1]-data[,i]
  }
  derta <- derta[-1,]
  derta[derta < 0] <- "-"
  derta[derta > 0] <- "+"
  i_number <- c(0)
  j = 1
  for (i in 1:nrow(derta)){
    A <- as.data.frame(derta[i,])
    S1 <- as.data.frame(seq[1:ncol(data)-1])
    S2 <- as.data.frame(seq[2:ncol(data)])
    A <- paste(A[1:nrow(A),],collapse = " ")
    S1 <- paste(S1[,1:ncol(S1)],collapse = " ")
    S2 <- paste(S2[,1:ncol(S2),],collapse = " ")
    if(A == S1 | A == S2){
      i_number[j] <- i
      j = j + 1
    }
  }
  data <- data[c(1,i_number+1),]
  data <- data.frame(t(data))
  return(data)
}##找出可能具有昼夜节律的通路（100%）
pmg <- datafilter(id1)
#write.csv(pmg,"h1.csv")
pmg <- datafilter(id2)
#write.csv(pmg,"h2.csv")
pmg <- datafilter(id3)
#write.csv(pmg,"h3.csv")
pmg <- datafilter(id4)
#write.csv(pmg,"h4.csv")
pmg <- datafilter(id5)
#write.csv(pmg,"h5.csv")

pmt <- read.csv("D:/Rfiles/rhythmst/pmt.csv", header = T, sep = ",")
la <- nrow(pmt)+1
n <- c(1)
j = 1
for (i in 1:nrow(pmt)){
  if(pmt[i,1] == 24){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
id1 <- data.frame(pmt[n[1]:(n[2]-1),])
id2 <- data.frame(pmt[n[2]:(n[3]-1),])
id3 <- data.frame(pmt[n[3]:(n[4]-1),])
id4 <- data.frame(pmt[n[4]:(n[5]-1),])
id5 <- data.frame(pmt[n[5]:(n[6]-1),])
sumzero <- function(x) sum(x == 0)
pmt <- datafilter(id1)
#write.csv(pmt,"h1.csv")
pmt <- datafilter(id2)
#write.csv(pmt,"h2.csv")
pmt <- datafilter(id3)
#write.csv(pmt,"h3.csv")
pmt <- datafilter(id4)
#write.csv(pmt,"h4.csv")
pmt <- datafilter(id5)
#write.csv(pmt,"h5.csv")

##可视化
par(mfrow=c(1,5))
plotlinec <- function(data,cols){
  ID <- rep(1:(ncol(data)-1),each = nrow(data))
  Time <- rep(data$Time, times = ncol(data)-1)
  lenr <- ncol(data)
  d <- c(0)
  for (i in 2:lenr){
    d <- c(d,data[,i])
  }
  Pp <- d[-1]
  pp <- data.frame(ID, Time, Pp)
  pp$ID<-as.numeric(pp$ID)
  pp$Time<-as.numeric(pp$Time)
  xrange<-range(pp$Time)
  yrange<-range(pp$Pp)
  plot(range(pp$Time),range(pp$Pp),type="n", xlab = "Sampling Time(h)", ylab = "Relative abundance(%)",xaxt="n")
  for (i in 1:(ncol(data)-1)){
    id <- subset(pp,ID == i)
    lines(id$Time,id$Pp,type = "b",lwd = 0.5,lty = 1,
          col = cols[i],pch = 19)
  }
  axis(side = 1, at = pp$Time, labels = pp$Time)
}
shared <- read.csv("D:/Rfiles/tetfig6new/shared.csv", header = T, sep = ",")
cols <- c("#006400","#96CCA8","#2E8B57","#C0EEDA","#7FFFAA","#447B66","#3CB371","#40E0D0","#CBD589","#86A845",
          "#32CD32","#B3C648","#FED46E","#FFDEAD","#F4A460","#DAA520","#FFFA32","#FFDF6F","#F4D160","#F0E68C",
          "#FFC86F","#8AC4D0","#5F9EA0","#58A3BC","#3E83A8","#3E64A8","#28527A","#134080","#323CC8","#1D6590",
          "#A8D4E0","#537496","#55A6CB","#A364A8","#A35AC8","#825AB4","#A378B5","#766092","#C71585","#DAB3DA",
          "#D986B1","#F2A6C2","#FBC1AD","#E164B4","#C71585","#EF978F",
          "#EFC2AD","#F56E4A","#F56E60","#F08080","#DB7093")
df <- data.frame(path = shared$path, col = cols,info = shared$info)##仅对shared Pathways上色
#h1--根据shared否进行整理后得到的数据
h1 <- read.csv("h1.csv", header = T, row.names = 1, sep = ",")
h1c <- read.csv("h1c.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h1c)
len2 <- ncol(h1)
c1 <- c(0)
for (i in 1:len1){
  c1[i] <-as.character(df[which(df$path==colnames(h1)[i]),2])
}
c1[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h1i <- read.csv("h1.csv", header = T, sep = ",")
plotlinec(h1i,c1)

#h2
h2 <- read.csv("h2c.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h2)
c2 <- c(0)
for (i in 1:len1){
  c2[i] <-as.character(df[which(df$path==colnames(h2)[i]),2])
}
h2i <- read.csv("h2c.csv", header = T, sep = ",")
plotlinec(h2i,c2)

#h3
h3 <- read.csv("h3g.csv", header = T, row.names = 1, sep = ",")
h3c <- read.csv("h3c.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h3c)
len2 <- ncol(h3)
c3 <- c(0)
for (i in 1:len1){
  c3[i] <-as.character(df[which(df$path==colnames(h3)[i]),2])
}
c3[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h3i <- read.csv("h3g.csv", header = T, sep = ",")
plotlinec(h3i,c3)

#h4
h4 <- read.csv("h4g.csv", header = T, row.names = 1, sep = ",")
h4c <- read.csv("h4c.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h4c)
len2 <- ncol(h4)
c4 <- c(0)
for (i in 1:len1){
  c4[i] <-as.character(df[which(df$path==colnames(h4)[i]),2])
}
c4[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h4i <- read.csv("h4g.csv", header = T, sep = ",")
plotlinec(h4i,c4)

#h5
h5 <- read.csv("h5g.csv", header = T, row.names = 1, sep = ",")
h5c <- read.csv("h5c.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h5c)
len2 <- ncol(h5)
c5 <- c(0)
for (i in 1:len1){
  c5[i] <-as.character(df[which(df$path==colnames(h5)[i]),2])
}
c5[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h5i <- read.csv("h5g.csv", header = T, sep = ",")
plotlinec(h5i,c5)

par(mfrow=c(1,5))
#h1rna
h1r <- read.csv("h1r.csv", header = T, row.names = 1, sep = ",")
h1cr <- read.csv("h1cr.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h1cr)
len2 <- ncol(h1r)
c1r <- c(0)
for (i in 1:len1){
  c1r[i] <-as.character(df[which(df$path==colnames(h1r)[i]),2])
}
c1r[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h1ir <- read.csv("h1r.csv", header = T, sep = ",")
plotlinec(h1ir,c1r)

#h2rna
h2r <- read.csv("h2r.csv", header = T, row.names = 1, sep = ",")
h2cr <- read.csv("h2cr.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h2cr)
len2 <- ncol(h2r)
c2r <- c(0)
for (i in 1:len1){
  c2r[i] <-as.character(df[which(df$path==colnames(h2cr)[i]),2])
}
c2r[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h2ir <- read.csv("h2r.csv", header = T, sep = ",")
plotlinec(h2ir,c2r)

#h3rna
h3r <- read.csv("h3r.csv", header = T, row.names = 1, sep = ",")
h3cr <- read.csv("h3cr.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h3cr)
len2 <- ncol(h3r)
c3r <- c(0)
for (i in 1:len1){
  c3r[i] <-as.character(df[which(df$path==colnames(h3cr)[i]),2])
}
c3r[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h3ir <- read.csv("h3r.csv", header = T, sep = ",")
plotlinec(h3ir,c3r)

#h4rna
h4r <- read.csv("h4r.csv", header = T, row.names = 1, sep = ",")
h4cr <- read.csv("h4cr.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h4cr)
len2 <- ncol(h4r)
c4r <- c(0)
for (i in 1:len1){
  c4r[i] <-as.character(df[which(df$path==colnames(h4cr)[i]),2])
}
c4r[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h4ir <- read.csv("h4r.csv", header = T, sep = ",")
plotlinec(h4ir,c4r)

#h5rna
h5r <- read.csv("h5r.csv", header = T, row.names = 1, sep = ",")
h5cr <- read.csv("h5cr.csv", header = T, row.names = 1, sep = ",")
len1 <- ncol(h5cr)
len2 <- ncol(h5r)
c5r <- c(0)
for (i in 1:len1){
  c5r[i] <-as.character(df[which(df$path==colnames(h5cr)[i]),2])
}
c5r[(len1+1):len2]<- rep("#D9D9D9", times=(len2-(len1+1)))
h5ir <- read.csv("h5r.csv", header = T, sep = ",")
plotlinec(h5ir,c5r)

###Species
##计算阶段
smg <- read.csv("D:/Rfiles/rhythmsg/smg.csv", header = T, sep = ",")
la <- nrow(smg)+1
n <- c(1)
j = 1
for (i in 1:nrow(smg)){
  if(smg[i,1] == 24){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
id1 <- data.frame(smg[n[1]:(n[2]-1),])
id2 <- data.frame(smg[n[2]:(n[3]-1),])
id3 <- data.frame(smg[n[3]:(n[4]-1),])
id4 <- data.frame(smg[n[4]:(n[5]-1),])
id5 <- data.frame(smg[n[5]:(n[6]-1),])
sumzero <- function(x) sum(x == 0)
smg <- datafilter(id1)
#write.csv(smg,"h1.csv")
smg <- datafilter(id2)
#write.csv(smg,"h2.csv")
smg <- datafilter(id3)
#write.csv(smg,"h3.csv")
smg <- datafilter(id4)
#write.csv(smg,"h4.csv")
smg <- datafilter(id5)
#write.csv(smg,"h5.csv")

smt <- read.csv("D:/Rfiles/rhythmst/smt.csv", header = T, sep = ",")
la <- nrow(smt)+1
n <- c(1)
j = 1
for (i in 1:nrow(smt)){
  if(smt[i,1] == 24){
    n[j] <- i
    j = j + 1
  }
}
n[j] <- la
id1 <- data.frame(smt[n[1]:(n[2]-1),])
id2 <- data.frame(smt[n[2]:(n[3]-1),])
id3 <- data.frame(smt[n[3]:(n[4]-1),])
id4 <- data.frame(smt[n[4]:(n[5]-1),])
id5 <- data.frame(smt[n[5]:(n[6]-1),])
sumzero <- function(x) sum(x == 0)
smt <- datafilter(id1)
#write.csv(smt,"h1.csv")
smt <- datafilter(id2)
#write.csv(smt,"h2.csv")
smt <- datafilter(id3)
#write.csv(smt,"h3.csv")
smt <- datafilter(id4)
#write.csv(smt,"h4.csv")
smt <- datafilter(id5)
#write.csv(smt,"h5.csv")

##可视化
par(mfrow = c(1,5))
plotlinec <- function(data,cols){
  ID <- rep(1:(ncol(data)-1),each = nrow(data))
  Time <- rep(data$Time, times = ncol(data)-1)
  lenr <- ncol(data)
  d <- c(0)
  for (i in 2:lenr){
    d <- c(d,data[,i])
  }
  Pp <- d[-1]
  pp <- data.frame(ID, Time, Pp)
  pp$ID<-as.numeric(pp$ID)
  pp$Time<-as.numeric(pp$Time)
  xrange<-range(pp$Time)
  yrange<-range(pp$Pp)
  plot(range(pp$Time),range(pp$Pp),type="n", xlab = "Sampling Time(h)", ylab = "Relative abundance(%)",xaxt="n")
  for (i in 1:(ncol(data)-1)){
    id <- subset(pp,ID == i)
    lines(id$Time,id$Pp,type = "b",lwd = 0.75,lty = 1,
          col = cols[i],pch = 19)
  }
  axis(side = 1, at = pp$Time, labels = pp$Time)
}##
species <- read.csv("D:/Rfiles/tetfig6new/species.csv", header = T, sep = ",")
cols <- c("#006400","#2E8B57","#C0EEDA","#447B66","#3CB371","#40E0D0","#CBD589","#86A845",
          "#32CD32","#B3C648","#FFDEAD","#F4A460","#DAA520","#FFDF6F","#F4D160","#F0E68C",
          "#FFC86F","#8AC4D0","#5F9EA0","#3E83A8","#3E64A8","#28527A","#134080","#323CC8","#1D6590",
          "#A8D4E0","#537496","#55A6CB","#A364A8","#825AB4","#A378B5","#766092","#C71585","#DAB3DA",
          "#D986B1","#F2A6C2","#FBC1AD","#E164B4","#C71585","#EF978F",
          "#EFC2AD","#F56E4A","#F56E60","#DB7093")
df <- data.frame(species = species$species, col = cols)#全部上色
#h1dna
h1s <- read.csv("h1s.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h1s)
c1s <- c(0)
for (i in 1:len){
  c1s[i] <-as.character(df[which(df$species==colnames(h1s)[i]),2])
}
h1si <- read.csv("h1s.csv", header = T, sep = ",")
plotlinec(h1si,c1s)

#h2dna
h2s <- read.csv("h2s.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h2s)
c2s <- c(0)
for (i in 1:len){
  c2s[i] <-as.character(df[which(df$species==colnames(h2s)[i]),2])
}
h2si <- read.csv("h2s.csv", header = T, sep = ",")
plotlinec(h2si,c2s)

#h3dna
h3s <- read.csv("h3s.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h3s)
c3s <- c(0)
for (i in 1:len){
  c3s[i] <-as.character(df[which(df$species==colnames(h3s)[i]),2])
}
h3si <- read.csv("h3s.csv", header = T, sep = ",")
plotlinec(h3si,c3s)

#h4dna
h4s <- read.csv("h4s.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h4s)
c4s <- c(0)
for (i in 1:len){
  c4s[i] <-as.character(df[which(df$species==colnames(h4s)[i]),2])
}
h4si <- read.csv("h4s.csv", header = T, sep = ",")
plotlinec(h4si,c4s)

#h5dna
h5s <- read.csv("h5s.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h5s)
c5s <- c(0)
for (i in 1:len){
  c5s[i] <-as.character(df[which(df$species==colnames(h5s)[i]),2])
}
h5si <- read.csv("h5s.csv", header = T, sep = ",")
plotlinec(h5si,c5s)

#h3rna
h3sr <- read.csv("h3sr.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h3sr)
c3sr <- c(0)
for (i in 1:len){
  c3sr[i] <-as.character(df[which(df$species==colnames(h3sr)[i]),2])
}
h3sri <- read.csv("h3sr.csv", header = T, sep = ",")
plotlinec(h3sri,c3sr)

#h4rna
h4sr <- read.csv("h4sr.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h4sr)
c4sr <- c(0)
for (i in 1:len){
  c4sr[i] <-as.character(df[which(df$species==colnames(h4sr)[i]),2])
}
h4sri <- read.csv("h4sr.csv", header = T, sep = ",")
plotlinec(h4sri,c4sr)

#h5rna
h5sr <- read.csv("h5sr.csv", header = T, row.names = 1, sep = ",")
len <- ncol(h5sr)
c5sr <- c(0)
for (i in 1:len){
  c5sr[i] <-as.character(df[which(df$species==colnames(h5sr)[i]),2])
}
h5sri <- read.csv("h5sr.csv", header = T, sep = ",")
plotlinec(h5sri,c5sr)

#画legend
par(mfrow=c(1,1))
plot(1,1)
legend('topright',col = cols, cex = .8, pch = 16,lty = 1,legend = df$info, ncol=2)
legend('topright',col = cols, cex = .8, pch = 16,lty = 1,legend = df$species, ncol=2)

