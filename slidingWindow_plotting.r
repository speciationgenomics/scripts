# Plot sliding window statistics 
#########################################################################

# Prepare input files:

# Read in the file with sliding window estimates of FST, pi and dxy
windowStats<-read.table("Pundamilia_kivu_div_stats.csv",header=T,sep=",")

# Remove windows with less than 10 kb covered
windowStats<-windowStats[windowStats$sites>10000,]


# Get chrom ends for making positions additive
chrom<-read.table("chrEnds.txt",header=T)
chrom$add<-c(0,cumsum(chrom$end)[1:21])
chrom$add<-chrom$add+c(0,cumsum(rep(10000000,times=21)))

# Make the positions of the divergence and diversity window estimates additive
windowStats$mid<-windowStats$mid+chrom[match(windowStats$scaffold,chrom$chr),3]


# Get fd estimates of 20 kb windows for NyerMak into NyerPyt
fd<-read.delim(file = "Pundamilia_ABBABABA.csv",sep=",",header=T,na.strings = "NaN")

# Remove the weird cases of negative fd or fd above 1
fd[(fd$fd<0 | fd$fd>1) & !is.na(fd$fd),"fd"]<-NA

# Make a density plot to show the distribution of number of sites used per window
# And to find a good threshold
plot(density(fd$sitesUsed,na.rm=T))
abline(v=10)

# Set windows with less than 5 sites to NA for all hybridization stats
fd[fd$sitesUsed<10&!is.na(fd$fd),c("fd","fdM","D")]<-c(NA,NA,NA)

# make positions additive
fd$mid<-fd$mid+chrom[match(fd$scaffold,chrom$chr),3]

##########################################################################
# Plotting:

# Combine 6 plots into a single figure:
par(mfrow=c(5,1),mar=c(0,5,0,1),cex.lab=1,cex.axis=1)

# Plot Fst between species at Makobe Island
plot(windowStats$mid,windowStats$Fst_NyerMak_PundMak,cex=0.5,pch=19,xaxt="n",
     ylab="Fst Makobe",ylim=c(0,1))
abline(h=mean(windowStats$Fst_NyerMak_PundMak,na.rm=T),col="grey")

# Plot Fst between species at Python Island
plot(windowStats$mid,windowStats$Fst_PundPyt_NyerPyt,cex=0.5,pch=19,xaxt="n",
     ylab="Fst Python",ylim=c(0,1))
abline(h=mean(windowStats$Fst_PundPyt_NyerPyt,na.rm=T),col="grey")

# Plot Dxy between species at Makobe Island
plot(windowStats$mid,windowStats$dxy_NyerMak_PundMak,cex=0.5,pch=19,xaxt="n",
     ylab="dxy Makobe",ylim=c(0,0.03))
abline(h=mean(windowStats$dxy_NyerMak_PundMak,na.rm=T),col="grey")

# Plot Dxy between species at Python Island
plot(windowStats$mid,windowStats$dxy_PundPyt_NyerPyt,cex=0.5,pch=19,xaxt="n",
     ylab="dxy Python",ylim=c(0,0.03))
abline(h=mean(windowStats$dxy_PundPyt_NyerPyt,na.rm=T),col="grey")

# Plot fd
plot(fd$mid,fd$fdM,cex=0.5,pch=19,xaxt="n",
     ylab="fdM",ylim=c(0,1))
abline(h=mean(fd$fdM,na.rm=T),col="grey")

# Add the LG names to the center of each LG
axis(1,at=chrom$add+chrom$end/2,tick = F,labels = 1:22)


##############################################################################

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


# Predict local recombination rates from a linkage map for each window
map<-read.table("linkage-map.bed",header=T)

# Remove positions that are too close together or without cM information:
map<-map[!is.na(map$cM),]
map<-map[c(TRUE,diff(map$POS)>20000),]

# Add predicted recombination rates:
addRecombinationRate <- function(dataset) {
  
  # Add a new column for the recombination rates
  dataset$rec<-rep(0,times=length(dataset$scaffold))
  
  # Loop through the 22 linkage groups (~chromosomes)
  for(chromosome in paste0("chr",1:22)){
    recChr<-map[map$CHROM==chromosome,]
    dataChr<-dataset[dataset$scaffold==chromosome,]
    tryCatch(
      d<-smooth.spline(recChr$POS,recChr$cM,spar=0.7),
      error=function(e) e)
    
    # get first derivative of smoothing curve = recombination rate
    # and predict recombination rates in cM/Mb at sliding window mid positions
    rec<-stats:::predict.smooth.spline(d,dataChr$mid,deriv=1)$y*1000000
    
    # set predicted recombination rates below 0 to 0
    rec[rec<0]<-0
    
    # Add the predicted recombination rates to the windowStats dataframe
    dataset[dataset$scaffold==chromosome,"rec"]<-rec
  }    
  
  return(dataset)
}
windowStats<-addRecombinationRate(windowStats)



##############################################################################
# TWISST
# R script to plot TWISST results by Simon:

simple.loess.predict <- function(x, y, span, weights = NULL, max = NULL, min = NULL){
  y.loess <- loess(y ~ x, span = span, weights = weights)
  y.predict <- predict(y.loess,x)
  if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
  if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
  y.predict
}

smooth_df <- function(x, df, span, col.names=NULL, weights=NULL, min=NULL, max=NULL){
  smoothed <- df
  if (is.null(col.names)){col.names=colnames(df)}
  for (col.name in col.names){
    print(paste("smoothing",col.name))
    smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, max = max, min = min, weights = weights)
  }
  smoothed
}

stack <- function(mat){
  upper <- t(apply(mat, 1, cumsum))
  lower <- upper - mat
  list(upper=upper,lower=lower)
}

interleave <- function(x1,x2){
  output <- vector(length= length(x1) + length(x2))
  output[seq(1,length(output),2)] <- x1
  output[seq(2,length(output),2)] <- x2
  output
}


plot_weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                         ylab="Weights", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
  #get x axis
  x = positions
  #if a two-column matrix is given - plot step-like weights with start and end of each window    
  if (is.matrix(x)==TRUE) {
    x = interleave(positions[,1],positions[,2])
    yreps=2
  }
  else {
    if (is.null(x)==FALSE) x = positions
    else x = 1:nrow(weights_dataframe)
    yreps=1
  }
  
  #set x limits
  if(is.null(xlim)) xlim = c(min(x), max(x))
  
  #if not adding to an old plot, make a new plot
  if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
  
  if (stacked == TRUE){
    y_stacked <- stack(weights_dataframe)
    for (n in 1:ncol(weights_dataframe)){
      y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
      y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
      polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], border=NA)
    }
  }
  else{
    for (n in 1:ncol(weights_dataframe)){
      y = rep(weights_dataframe[,n],each=yreps)
      polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA)
      lines(x,y, type = "l", col = line_cols[n])
    }
  }
}

options(scipen = 5)

########### input data ################

# Input data for piscivore vs detritivores in Kagera and Victoria
weights_file <- "D:/Dropbox/victoriaGenomes/TWISST/piscDetrViKa.weights.csv.gz"
window_data_file <- "D:/Dropbox/victoriaGenomes/TWISST/piscDetrViKa.data.tsv"

########## read data ##################
weights = read.table(gzfile(weights_file), header = F)  # if gzipped
names(weights)<-c("topo1","topo2","topo3")

#normalise rows so weights sum to 1
weights <- weights / apply(weights, 1, sum)
#retrieve the names of the topologies
topoNames = names(weights)

window_data = read.table(window_data_file, header = T)

#exclude any rows where data is missing
good_rows = which(is.na(apply(weights,1,sum)) == F)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]

#exclude rows with too much missing data
good_rows = which(window_data$sites>100)
weights <- weights[good_rows,]
window_data = window_data[good_rows,]


# Make positions additive
window_data<-cbind(window_data,"nonAddStart"=window_data$start,"nonAddEnd"=window_data$end)
chrom<-read.csv("D:/Dropbox/MwanzaSpeciationTransect/whole_genome_sequencing/chromAddInfo.csv",sep=";")
for(chr in levels(droplevels(window_data$scaffold))){
  add<-chrom[chrom$chr==chr,"add"]
  window_data[window_data$scaffold==chr,"start"]<-window_data[window_data$scaffold==chr,"start"]+add
  window_data[window_data$scaffold==chr,"mid"]<-window_data[window_data$scaffold==chr,"mid"]+add
  window_data[window_data$scaffold==chr,"end"]<-window_data[window_data$scaffold==chr,"end"]+add
}

weights<-cbind(pos=window_data[,4],weights)


########### choose colours for plot ########
#some nice contrasted colours.
cols = c(
  "#993F00", #Uranium
  "#0000DC", #Blue
  "#808080", #grey
  "#993F00", #Caramel
  "#4C005C", #Damson
  "#191919", #Ebony
  "#F0A3FF", #Amethyst
  "#2BCE48", #Green
  "#FFCC99", #Honeydew
  "#808080", #Iron
  "#94FFB5", #Jade
  "#8F7C00", #Khaki
  "#9DCC00", #Lime
  "#C20088", #Mallow
  "#003380", #Navy
  "#FFA405", #Orpiment
  "#FFA8BB", #Pink
  "#426600", #Quagmire
  "#FF0010", #Red
  "#5EF1F2", #Sky
  "#00998F", #Turquoise
  "#740AFF", #Violet
  "#990000", #Wine
  "#FFFF80", #Xanthin
  "#FFFF00", #Yellow
  "#FF5005" #Zinnia
)






