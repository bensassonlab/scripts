
rm(list=ls())
pdf("P34048.nearest.pdf")
data<-read.table("P34048.newnearest.tsv",header=T)
attach(data)
head(data)
tail(data)

par(mfrow=c(8,1),mar=c(1,1,1,1),yaxt='n',xlab='',bty="n")



# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr1.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr1.mfa"]-100000,0,winpos[infile=="chr1.mfa"],100,col=as.vector(cladecolor[infile=="chr1.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr2.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr2.mfa"]-100000,0,winpos[infile=="chr2.mfa"],100,col=as.vector(cladecolor[infile=="chr2.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr3.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr3.mfa"]-100000,0,winpos[infile=="chr3.mfa"],100,col=as.vector(cladecolor[infile=="chr3.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr4.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr4.mfa"]-100000,0,winpos[infile=="chr4.mfa"],100,col=as.vector(cladecolor[infile=="chr4.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr5.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr5.mfa"]-100000,0,winpos[infile=="chr5.mfa"],100,col=as.vector(cladecolor[infile=="chr5.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr6.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr6.mfa"]-100000,0,winpos[infile=="chr6.mfa"],100,col=as.vector(cladecolor[infile=="chr6.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chr7.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chr7.mfa"]-100000,0,winpos[infile=="chr7.mfa"],100,col=as.vector(cladecolor[infile=="chr7.mfa"]))


# PRINT STRAIGHTFORWARD GENOME PLOTS USING USER-SPECIFIED COLOUR FOR BINS WITH TIED MATCHES
# data from P34048.newnearest.tsv
plot(c(0,max(winpos)),c(0,100),type="n",main="P34048.chrR.mfa",yaxt='n',xaxt='n')
rect(winpos[infile=="chrR.mfa"]-100000,0,winpos[infile=="chrR.mfa"],100,col=as.vector(cladecolor[infile=="chrR.mfa"]))


# OVERWRITE WITH SMALLER RECTANGLES TO SHOW MULTIPLE COLORS IN CASE OF TIES 
# data is added in using perl and is not contained in P34048.newnearest.tsv
rect(400000-100000,0,400000,50,col="blue")


# OVERWRITE WITH SMALLER RECTANGLES TO SHOW MULTIPLE COLORS IN CASE OF TIES 
# data is added in using perl and is not contained in P34048.newnearest.tsv
rect(400000-100000,50,400000,100,col="#5aae61")
