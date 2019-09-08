#!/usr/bin/Rscript

# (c) Joana Meier 2018
# R-script to plot the model
# requires <prefix>.bestlhoods, <prefix>.tpl <prefix>_maxL.par 
# ideally also a <prefix>.AIC file

# Usage: Rscript plotModel.r -p <prefix> -l <pop0,pop1>

# Load libaries
library(optparse)
library(shape)
library(plotrix)

# Read input arguments
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="infile prefix", metavar="character"),
  make_option(c("-l", "--list"), type="character", default=NULL,
              help="list of populations", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check whether input prefix is provided
# if not provided, abort
if(is.null(opt$prefix)){
  stop(paste("Aborted. Please provide a prefix with -p <prefix>.\nUsage: plotModel.r -p prefix",sep=""))
}
if(is.null(opt$list)){
  stop("Aborted. Please provide a comma-separated list of population names e.g. -l pop1,pop2")
}

prefix=opt$prefix
populations=unlist(strsplit(opt$list,","))
popN=length(populations)
path=""

# requires R packages shape and plotrix
# Function to plot the model based on *maxL.par
plotModelBestParams <- function (prefix="",logNe=FALSE) {
  # read the input files (tpl and bestlhoods files)
  tpl<-readLines(paste(path,prefix,".tpl",sep=""),skip=3,)
  maxLpar<-readLines(paste(path,prefix,"_maxL.par",sep=""))
  
  # parse the tpl file
  NeVar<-tpl[(grep("Population effective sizes",tpl)+1):(grep("//Haploid",tpl)-1)]
  Ne<-as.integer(maxLpar[(grep("Population effective sizes",maxLpar)+1):(grep("//Haploid",maxLpar)-1)])
  orNe<-Ne
  totN<-sum(Ne)
  eventsVar<-tpl[(grep("//historical event",tpl)+2):(grep("//Number of independent loci",tpl)-1)]
  events<-maxLpar[(grep("//historical event",maxLpar)+2):(grep("//Number of independent loci",maxLpar)-1)]
  splitsVar<-strsplit(x=eventsVar,split=" ")
  events<-sub(pattern = "  ",replacement = " ",x = events)
  splits<-strsplit(x=events,split=" ")
  age<-rep(0,times=popN)
  migMatrixN<-tpl[(grep("//Number of migration matrices",tpl)+1)]
  maxMig<-0
  
  # remove events that do not merge populations (e.g. change of migration matrix) (same dummy population is source and sink)
  splits[sapply(splits, "[[", 2)==sapply(splits, "[[", 3)]<-NULL
  splitsVar[sapply(splitsVar, "[[", 2)==sapply(splitsVar, "[[", 3)]<-NULL
  
  
  # if >0 matrices in the tpl file, read the first one
  edgeN=1
  mig<-NA
  migTimes<-NA
  for(i in 0:(as.integer(migMatrixN)-1)){
    lN<-grep(paste("//Migration matrix ",i,sep=""),maxLpar)
    migMatrix<-maxLpar[(lN+1):(lN+popN)]
    if(i<1) mig<-strsplit(x=migMatrix,split="\\s+")
    if(i>0) mig<-c(mig,strsplit(x=migMatrix,split="\\s+"))
  }
  migProp<-as.double(unlist(mig)[unlist(mig)!="0"])
  edgeN<-length(migProp)
  maxMig<-max(migProp)
  
  
  # add real split times and sort by split times
  for(i in 1:length(splits)){
    splitName<-as.character(splitsVar[[i]][1])
    splits[[i]]<-c(splitName,splits[[i]])
  }
  splits<-splits[order(as.integer(sapply(splits,head,2)[2,]),decreasing=F)]
  
  # set age of populations (up to Ne change)
  maxT<-0
  maxN<-max(Ne)
  adjNe<-Ne
  for(i in 1:length(splits)){
    splitName<-splits[[i]][1]
    time<-as.integer(splits[[i]][2])
    if(time>maxT) maxT<-time
    sourcePop<-as.integer(splits[[i]][3])+1
    propMerge<-as.integer(splits[[i]][5])
    AncN<-as.double(splits[[i]][6])
    sinkPop<-as.integer(splits[[i]][4])+1
    if(propMerge==1 && age[sourcePop]==0) age[sourcePop]=time
    if(AncN!=1){
      if(age[sinkPop]==0) age[sinkPop]=time
      # to get maxN for plotting, account for ancestral Ne changes
      newSize<-adjNe[sinkPop]*AncN
      adjNe[sinkPop]<-newSize
      if(newSize>maxN) maxN=newSize
    }
  }
  
  # find oldest time point to set plotting limits
  maxT<-as.integer(maxT*1.2)
  firstSplit<-min(age)
  age[age==0]<-maxT
  migSpace<-maxT/20*edgeN
  div<-maxN*2
  
  # get the plotting area
  par(mfrow=c(1,1),mai=c(1,1,1,0.5),mgp=c(1.5,0.2,0))
  plot(x=c(1,popN),y=0:1, xlim=c(-0.1,(popN+0.5)),ylim=c(-migSpace,maxT),
       type = "n",xlab="",ylab="time (generations)",xaxt="n",yaxt="n",
       main=prefix)
  axis(2,at=pretty(0:maxT))
  
  # add population bodies
  for(i in 1:popN){
    if(logNe){
      div<-log(maxN)*2
      rect(xleft=i-log(Ne[i])/div,xright=i+log(Ne[i])/div,
           ybottom=0,ytop=age[i],col="grey",border=NA)
    }
    else{
      if(Ne[i]/div>0.01) rect(xleft=i-Ne[i]/div,xright=i+Ne[i]/div,
                              ybottom=0,ytop=age[i],col="grey",border=NA)
      else rect(xleft=i-0.01,xright=i+0.01,
                ybottom=0,ytop=age[i],col="grey",border=NA)
    }
  }
  text(x=c(1:popN),y=0,labels=NeVar,adj=c(0.5,0))
  text(x=c(1:popN),y=0,labels=Ne,adj=c(0.5,1))
  
  
  # if population sizes change at split, add rectangles
  for(i in 1:length(splits)){
    ancResize<-as.double(splits[[i]][6])
    timeResize<-as.integer(splits[[i]][2])
    
    # if change of ancestral Ne, draw new rectangle
    if(ancResize!=1){
      popPos=as.integer(splits[[i]][4])+1
      newN=as.double(Ne[popPos])*ancResize
      # draw a white rectangle hiding potential Ne rectangles
      if(logNe){
        rect(xleft=popPos-log(Ne[popPos])/div,
             xright=popPos+log(Ne[popPos])/div,lwd=1,
             ybottom=timeResize,ytop=maxT,col="white",border="white")
      }else{
        if(Ne[popPos]/div<0.01){rect(xleft=popPos-0.01,xright=popPos+0.01,
                                     ybottom=timeResize,ytop=maxT,col="white",border="white",lwd=2)
        }else rect(xleft=popPos-(Ne[popPos]/div),xright=popPos+(Ne[popPos]/div),
                   ybottom=timeResize,ytop=maxT,col="white",border="white",lwd=1)
      }
      # change Ne to the new value
      Ne[popPos]<-newN
      
      # draw the rectangle with the correct new Ne
      if(logNe){
        rect(xleft=popPos-log(newN)/div,
             xright=popPos+log(newN)/div,ybottom=timeResize,
             ytop=maxT,col="grey",border="grey")
      }else{
        if(newN/div<0.01){
          rect(xleft=popPos-0.01,xright=popPos+0.01,
               ybottom=timeResize,ytop=maxT,col="grey",border=NA)
        }else rect(xleft=popPos-newN/div,xright=popPos+newN/div,
                   ybottom=timeResize,ytop=maxT,col="grey",border=NA)
      }
      # add text
      txt<-round(ancResize,digits=2)
      #text(x=popPos,y=timeResize,labels=paste(txt,"x",sep=""),adj=c(0.5,1))
      text(x=popPos,y=timeResize,labels=round(Ne[popPos],digits=0),adj=c(0.5,0))
    }
    
    # reduce age of source pop if merging into other pop
    if(as.integer(splits[[i]][5])==1){
      sourcePop<-as.integer(splits[[i]][3])+1
      if(logNe) rect(xleft=sourcePop-(log(Ne[sourcePop])/div),
                     xright=sourcePop+(log(Ne[sourcePop])/div),
                     ybottom=timeResize,ytop=maxT,col="white",border="white",lwd=1)
      else{
        if(Ne[sourcePop]/div<0.01){ rect(xleft=sourcePop-0.01,xright=sourcePop+0.01,
                                         ybottom=timeResize,ytop=maxT,col="white",border="white",lwd=1)
        }else rect(xleft=sourcePop-(Ne[sourcePop]/div),xright=sourcePop+(Ne[sourcePop]/div),
                   ybottom=timeResize,ytop=maxT,col="white",border="white",lwd=1)
      }
    }
    splits[[i]]
  }
  
  # draw arrows showing which populations merge and add split times to y axis
  splitTimes<-matrix(nrow=length(splits),ncol=2)
  require("shape")
  for(i in 1:length(splits)){
    splitName<-splits[[i]][1]
    time<-as.integer(splits[[i]][2])
    sourcePop<-as.integer(splits[[i]][3])+1
    sinkPop<-as.integer(splits[[i]][4])+1
    propMerge<-as.double(splits[[i]][5])
    Arrows(x0=sourcePop,x1=sinkPop,y0=time,y1=time,lwd=propMerge,
           arr.width=0.1,arr.adj=1,arr.type="simple")
    if(propMerge<1) text(x=mean(c(sourcePop,sinkPop)),y=time*1.2,
                         labels=round(propMerge,digits = 3))
    splitTimes[i,]<-c(time,paste(splitName,time))
    points(x=sourcePop,y=time,pch=20)
  }
  axis(side=2,hadj=-0.2,at=as.integer(splitTimes[,1]),
       labels=splitTimes[,2],las=1,tcl=0.5)
  
  # add a title for the split times
  text(x=-0.25,y=maxT,labels="Split times",adj=c(0,0))
  
  # add migration edges (in first migration matrix)
  if(migMatrixN>0){
    
    # y position for first arrow
    add=-migSpace
    
    maxMig=0.015
    
    
    # plot migration edges
    require(plotrix)
    for(line in 1:length(mig)){
      edges<-as.double(mig[[line]])
      for(sinkPop in 1:length(edges)){
        edge<-edges[sinkPop]
        if(edge!="0"){
          standStrength<-edge/maxMig
          sourcePop<-line%%popN
          if(sourcePop==0) sourcePop=popN
          Arrows(x0=sourcePop,x1=sinkPop,y0=add,y1=add,lwd=1,
                 arr.width=0.1,arr.adj=1,arr.type="simple",
                 col=rgb(1,1-standStrength,0,maxColorValue=1))
          text(x=sourcePop,y=-0.2+add,cex=0.8,
               labels=paste(format(edge,scientific=TRUE,digits=3),"/",
                            round(edge*orNe[sinkPop],digits = 3)))
          add=add+((migSpace*0.8)/edgeN)
        }
      }
    }
    lut=colorRampPalette(c("yellow","red"))(50)
    nticks=edgeN/2
    scale = (length(lut)-1)/(migSpace)
    ticks<-format(seq(0,maxMig,len=nticks),scientific=T,digits=2)
    axis(2,at=seq(-migSpace*1.2,-migSpace*0.2,len=length(ticks)),
         labels=ticks,las=1,cex=0.8,tcl=0.3)
    for (i in 1:(length(lut)-1)) {
      y = (i-1)/scale-(migSpace*1.2)
      rect(-0.25,y,0.1,y+1/scale, col=lut[i], border=NA)
    }
    text(x=-0.1,y=-migSpace*0.2,adj=c(0.5,-0.2),labels="mig")
  }
}

pdf(paste(path,prefix,"_model.pdf",sep=""),paper="a4",height=20,width=10)

# plot the model with the bestlhoods parameters
par(mfrow=c(1,1),oma=c(10,0,10,0),cex.axis=0.8,cex.main=1,cex.lab=1,cex=0.8)
plotModelBestParams(prefix,logNe=F)

# add info about likelihood
if(file.exists(paste(path,prefix,"/",prefix,".bestlhoods",sep=""))){
  bestlk<-read.table(paste(path,prefix,"/",prefix,".bestlhoods",sep=""),header=T,sep="\t")
}else{
  bestlk<-read.table(paste(path,prefix,".bestlhoods",sep=""),header=T,sep="\t")
}
if(file.exists(paste0(path,prefix,"/",prefix,".AIC"))){
  aic<-read.table(paste0(path,prefix,"/",prefix,".AIC"),header=T,sep="\t")
}else if(file.exists(paste0(path,prefix,".AIC"))){
  aic<-read.table(paste0(path,prefix,".AIC"),header=T,sep="\t")
}else{
  aic<-NA
}
title(sub=paste("MaxEstLhood: ",round(bestlk["MaxEstLhood"],digits=1),
                ", MaxObsLhood: ",round(bestlk["MaxObsLhood"],digits=1),
                ", diff: ",round(aic$deltaL,digits=1),
                ", AIC: ",round(aic$AIC,digits=1),"\n","Note: All parameters (incl. migration) are plotted backward in time and in haploid numbers",sep=""))

dev.off()
