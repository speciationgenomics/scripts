#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r <prefix> <infofile>

# Read arguments from the command line
args=commandArgs(TRUE)

# Assign the first argument to prefix
prefix=args[1]

# Get individual names in the correct order
labels<-read.table(args[2],sep=" ")

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
labels$n<-as.factor(labels$pop)
levels(labels$n)<-c(4,2,5,1,3)
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
tbl2=read.table(paste0(prefix,".2.Q"))
tbl3=read.table(paste0(prefix,".3.Q"))
tbl4=read.table(paste0(prefix,".4.Q"))
tbl5=read.table(paste0(prefix,".5.Q"))

# Prepare spaces between the species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.5)}
spaces<-spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual
tiff(file=paste0(prefix,".tiff"),width = 1000, height = 600,res=300)
par(mfrow=c(4,1),mar=c(0,1,0,0),oma=c(2,1,9,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
bp<-barplot(t(as.matrix(tbl2[order(labels$n),])), col=rainbow(n=2),xaxt="n", border=NA,ylab="K=2",yaxt="n",space=spaces)
axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)
barplot(t(as.matrix(tbl3[order(labels$n),])), col=rainbow(n=3),xaxt="n", border=NA,ylab="K=3",yaxt="n",space=spaces)
barplot(t(as.matrix(tbl4[order(labels$n),])), col=rainbow(n=4),xaxt="n",  border=NA,ylab="K=4",yaxt="n",space=spaces)
barplot(t(as.matrix(tbl5[order(labels$n),])), col=rainbow(n=5),xaxt="n", xlab="Individual #", border=NA,ylab="K=5",yaxt="n",space=spaces)
dev.off()
