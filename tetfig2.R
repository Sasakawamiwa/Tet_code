setwd("D:/Rfiles/tetfig2")
library(vegan)
library(ellipse)

sp <- data.frame(t(read.csv("spceispcoa.csv",header = T,row.names = 1,sep = ",")))
pp <- data.frame(t(read.csv("pathpcoa.csv",header = T,row.names = 1,sep = ",")))
g <- read.csv("pcoagroup.csv", header = T,sep = ",")
g$group <- factor(g$group, levels=c("PreFD","PreFR","PostFD","PostFR","PostID","PostIR"))
group <- g$group

###传统画法
pcoaplot2 <- function(data,g){
  group <- g$group
  bray <- vegdist(data,"bray")
  jaccard <- vegdist(data,"jaccard")
  data_cca <- cca(data)
  data_cap <- capscale(bray ~ -1)
  pcoa<- capscale(bray ~ group)
  xl = pcoa$CA$eig[1]/pcoa$tot.chi*100
  yl = pcoa$CA$eig[2]/pcoa$tot.chi*100
  xlab = paste("PCoA axis1:", round(xl[1], 2), "%"," Variation", sep = "")
  ylab = paste("PCoA axis2:", round(yl[1], 2), "%"," Variation", sep = "")
  sum <- summary(pcoa)
  site <- sum$sites[,1:2]
  site <- data.frame(site)
  colnames(site) <- c('PCoA1', 'PCoA2')
  site$group <- g$group
  cols <- c("#FFD700","#FFDF6F","#FBC1AD","#D986B1","#7DABD0","#1D6590")
  plot(site$PCoA1, site$PCoA2, type = "n", xlab = xlab, ylab = ylab)
  points(site$PCoA1, site$PCoA2, col = cols[g$group], cex = 0.8, pch = 20)
  ordiellipse(pcoa,group,lty=2,col=cols)
  legend('topright',pch=20,col=cols,legend =levels(group))
}
par(mfrow=c(1,2))
pcoaplot2(sp,g)
pcoaplot2(pp,g)
