rm(list = ls())

# Prepare input files:
# Read in the file with sliding window estimates of FST, pi and dxy
windowStats<-read.csv("./genome_scan/Pundamilia_kivu_div_stats.csv",header=T)

# Let's have a look at what columns are present
names(windowStats)
length(windowStats$scaffold)

# Get chrom ends for making positions additive
chrom<-read.table("./genome_scan/chrEnds.txt",header=T)
chrom$add<-c(0,cumsum(chrom$end)[1:21])

# Make the positions of the divergence and diversity window estimates additive
windowStats$mid<-windowStats$mid+chrom[match(windowStats$scaffold,chrom$chr),3]

# Get fd estimates of 20 kb windows for NyerMak into NyerPyt
fd<-read.delim(file = "./genome_scan/Pundamilia_ABBABABA.csv",sep=",",header=T,na.strings = "NaN")

# Remove windows without mid position (windows with 0 sites)
fd<-fd[!is.na(fd$mid),]

# make positions additive
fd$mid<-fd$mid+chrom[match(fd$scaffold,chrom$chr),3]

##########################################################################
# Plotting statistics along the genome:

# Combine 2 plots into a single figure:
par(mfrow=c(2,1),oma=c(3,0,0,0),mar=c(1,5,1,1))

# Plot Fst between species at Makobe Island
plot(windowStats$mid,windowStats$Fst_NyerMak_PundMak,cex=0.5,pch=19,xaxt="n",
     ylab="Fst Makobe",ylim=c(0,1),col=windowStats$scaffold)
abline(h=mean(windowStats$Fst_NyerMak_PundMak,na.rm=T),col="grey")

# Plot Fst between species at Python Island
plot(windowStats$mid,windowStats$Fst_PundPyt_NyerPyt,cex=0.5,pch=19,xaxt="n",
     ylab="Fst Python",ylim=c(0,1),col=windowStats$scaffold)
abline(h=mean(windowStats$Fst_PundPyt_NyerPyt,na.rm=T),col="grey")
# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = F,labels = 1:22)

# Plot Dxy between species at Makobe Island
plot(windowStats$mid,windowStats$dxy_NyerMak_PundMak,cex=0.5,pch=19,xaxt="n",
     ylab="dxy Makobe",ylim=c(0,0.03),col=windowStats$scaffold)
abline(h=mean(windowStats$dxy_NyerMak_PundMak,na.rm=T),col="grey")

# Plot Dxy between species at Python Island
plot(windowStats$mid,windowStats$dxy_PundPyt_NyerPyt,cex=0.5,pch=19,xaxt="n",ylab="dxy Python",ylim=c(0,0.03),col=windowStats$scaffold)
abline(h=mean(windowStats$dxy_PundPyt_NyerPyt,na.rm=T),col="grey")

# Plot fd along the genome
plot(fd$mid,fd$fdM,cex=0.5,pch=19,xaxt="n",
     ylab="fdM",ylim=c(0,1),col=windowStats$scaffold)
abline(h=mean(fd$fdM,na.rm=T),col="grey")

# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = F,labels = 1:22)

# Are dxy and fst correlated?
par(mfrow=c(1,1),mar=c(5,5,1,1),oma=c(1,1,1,1))

plot(windowStats$dxy_NyerMak_PundMak,windowStats$Fst_NyerMak_PundMak,
     xlab="dxy",ylab="Fst",pch=19,cex=0.3)
# ignoring the outliers:
plot(windowStats$dxy_NyerMak_PundMak,windowStats$Fst_NyerMak_PundMak,
     xlab="dxy",ylab="Fst",pch=19,cex=0.3,xlim=c(0,0.01))

# Is dxy correlated with pi, indicating that it is mostly driven by variance in mutation and recombination rates?
plot(rowMeans(cbind(windowStats$pi_PundMak,windowStats$pi_NyerMak),na.rm=T),
     windowStats$dxy_NyerMak_PundMak,cex=0.3,
     xlab="pi",ylab="dxy")
abline(a=0,b=1,col="grey")

# If we correct for mutation rate differences with dxy to kivu
plot(windowStats$dxy_NyerMak_PundMak/rowMeans(cbind(windowStats$dxy_NyerMak_kivu,
                                                    windowStats$dxy_PundMak_kivu),na.rm=T),
     windowStats$Fst_NyerMak_PundMak,ylab="Fst",xlab="dxy normalized",pch=19,cex=0.3)


##############################################################################

# Supplementary


# What are the FST distributions between the two species pairs?
par(mfrow=c(2,2))
boxplot(windowStats$Fst_PundPyt_NyerPyt,windowStats$Fst_NyerMak_PundMak)
abline(h=0)

# Are shared high Fst windows in regions of low recombination?
sharedHighFst<-windowStats[windowStats$Fst_NyerMak_PundMak>0.1 &
                             windowStats$Fst_PundPyt_NyerPyt>0.1, ]
normalFST<-windowStats[windowStats$Fst_NyerMak_PundMak<0.1 &
                             windowStats$Fst_PundPyt_NyerPyt<0.1,]

# Add the ratio of dxy between Makobe to dxy to kivu
sharedHighFst$ratio<-sharedHighFst$dxy_NyerMak_PundMak/sharedHighFst$dxy_NyerMak_kivu
sharedHighFst<-normalFST[!is.infinite(sharedHighFst$ratio)|is.na(sharedHighFst$ratio),]
normalFST$ratio<-normalFST$dxy_NyerMak_PundMak/normalFST$dxy_NyerMak_kivu
normalFST<-normalFST[!is.infinite(normalFST$ratio)|is.na(normalFST$ratio),]

# How many sites of shared high FST (>0.1) are there?
length(sharedHighFst$scaffold)

# Compare stats at windows with high fst to windows with Fst<0.1
vioplot(sharedHighFst$dxy_PundPyt_NyerPyt,normalFST$dxy_NyerMak_PundMak)
vioplot(sharedHighFst$pi_PundPyt,normalFST$pi_PundPyt)
vioplot(sharedHighFst$pi_NyerPyt,normalFST$pi_NyerPyt)
vioplot(sharedHighFst$pi_PundMak,normalFST$pi_PundMak)
vioplot(sharedHighFst$pi_NyerMak,normalFST$pi_NyerMak)

# Do shared high FST windows contain haplotypes predating the Victoria-Kivu split?
vioplot(sharedHighFst$ratio,normalFST$ratio,horizontal = T)

install.packages("beanplot")
require("beanplot")

par(mfrow=c(1,1),mar=c(1,3,1,1),mgp=c(1.6,0.5,0),cex.axis=0.9,cex=1.3,xaxs="i",yaxs="i")
bp<-beanplot(normalFST$ratio,sharedHighFst$ratio,log = "",
             side='both', border='NA',
             col=list('cornflowerblue','blue'),
             ylab='maximum relative divergence' ,what=c(0,1,0,0),xaxt="n",
             ylim=c(0,3),cex=2)
abline(h=1)
par(xpd=F)

# Function to plot the density of dots in a grid
colPlot <- function(dataset=stats,varx,vary,minx=0,miny=0,maxx=0.02,maxy=maxx,title="",xlab=varx,ylab=vary,corr=T){
  rcbpal<-c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026")
  data<-dataset
  rf<-colorRampPalette((rcbpal))
  r<-rf(12)
  xbin<-seq(minx,maxx,length.out = 100)
  ybin<-seq(miny,maxy,length.out = 100)
  # colx=grep(names(data),pattern=varx)
  # coly=grep(names(data),pattern=vary)
  colx=varx
  coly=vary
  freq<-as.data.frame(table(findInterval(x = data[,colx],vec = xbin,all.inside=T),findInterval(x = data[,coly],vec = ybin,all.inside = T)))
  freq[,1] <- as.numeric(as.character(freq[,1]))
  freq[,2] <- as.numeric(as.character(freq[,2]))
  freq2D<-matrix(0,nrow=length(xbin),ncol=length(ybin))
  freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
  image(xbin,ybin,log10(freq2D),col=r,breaks = seq(minx,max(log10(freq2D)),length.out=length(r)+1),
        xlab=xlab,ylab=ylab,xlim=c(minx,maxx),ylim=c(miny,maxy))
  barval<-round(10^(seq(minx,max(log10(freq2D)),length.out = length(r)+1)))
  title(main=title)
  lm<-summary(lm(data[,coly]~data[,colx]))
  r2<-round(lm$r.squared,3)
  if(corr) legend("top",legend=bquote(r^2 == .(r2)),bty = "n",cex=1.2)
  abline(lm(data[,coly]~data[,colx]))
}

par(mfrow=c(2,2),mar=c(5,5,0,0))

# Are the FST values of the two species pairs correlated?
colPlot(dataset = windowStats,varx = "Fst_PundPyt_NyerPyt",vary = "Fst_NyerMak_PundMak",maxx = 1,maxy=1)
colPlot(dataset = windowStats,varx = "dxy_PundPyt_NyerPyt",vary = "dxy_NyerMak_PundMak",maxx = 0.03,maxy=0.03)

# Are dxy and pi correlated?
colPlot(dataset = windowStats,varx = "pi_NyerPyt",vary = "dxy_PundPyt_NyerPyt",maxx = 0.03,maxy=0.03)
colPlot(dataset = windowStats,varx = "pi_NyerMak",vary = "dxy_NyerMak_PundMak",maxx = 0.03,maxy=0.03)

