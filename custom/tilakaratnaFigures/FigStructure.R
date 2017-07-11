# Edited Rcmds output from structurePrint.pl

pdf('figStructure.pdf')
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

data<-read.table("allCensFrom_gfa2.infile.k6.structure.data")
attach(data)
x<-rbind(data[,6],data[,7],data[,8],data[,9],data[,10],data[,11])
par(mar=c(7, 4, 4, 2) + 0.1,cex=0.6,xpd=T)
#z<-barplot(x,beside=FALSE,col=c("#8c510a","#4daf4a","#377eb8","#e41a1c","orange","#ffff33"),names.arg=data[,2],las=2, border=NA,space=0,main='a) Worldwide, K=6',xaxs = "i")
z<-barplot(x,beside=FALSE,col=c("#8c510a","#4daf4a","#377eb8","#e41a1c","orange","#ffff33"),names.arg=data[,2],las=2, border=NA,space=0,main='a) Worldwide strains',xaxs = "i")

# show population names and numbers
par(cex=1)
text(10,0.5,"Europe")
text(50,0.5,"USA")
text(66.5,0.5,"1")
text(69,0.5,"2")
text(78,0.5,"3")

# show ref strains
text(z[which(data[,2]=="UWOPS03_461")],-0.04,"*")
text(z[which(data[,2]=="UWOPS05_217")],-0.04,"*")
#text(z[which(data[,2]=="UWOPS05_227")],-0.04,"*")
#text(z[which(data[,2]=="Y12m")],-0.04,"*")
text(z[which(data[,2]=="Y9m")],-0.04,"*")
text(z[which(data[,2]=="K11m")],-0.04,"*")
text(z[which(data[,2]=="NCYC110m")],-0.04,"*")
text(z[which(data[,2]=="DBVPG6044m")],-0.04,"*")

text(z[which(data[,2]=="L1374m")],-0.04,"*")
text(z[which(data[,2]=="DBVPG1106m")],-0.04,"*")
text(z[which(data[,2]=="ZP562s1")],-0.04,"*")
text(z[which(data[,2]=="ZP560s1")],-0.04,"*")
text(z[which(data[,2]=="YPS128m")],-0.04,"*")
text(z[which(data[,2]=="YPS606m")],-0.04,"*")
text(z[which(data[,2]=="SDO1s1")],-0.04,"*")
text(z[which(data[,2]=="SDO3s1")],-0.04,"*")


# Figure B
data<-read.table("allCensFrom_gfa2.infile.k3.structure.data")
attach(data)
x<-rbind(data[,6],data[,7],data[,8])
par(mar=c(7, 7, 4, 5) + 0.1,cex=0.6)
#barplot(x,beside=FALSE,col=c("#762a83","#8c510a","#a6dba0"),names.arg=data[,2],las=2, border=NA,space=0,main='b) Wine or European oaks, K=3',xaxs = "i")
barplot(x,beside=FALSE,col=c("#762a83","#8c510a","#a6dba0"),names.arg=data[,2],las=2, border=NA,space=0,main='b) Wine or European oak strains',xaxs = "i")


# show population names and numbers
par(cex=1)
text(7,0.5,"4")
text(18.5,0.5,"5")

# show ref strains
text(z[which(data[,2]=="L1374m")],-0.04,"*")
text(z[which(data[,2]=="DBVPG1106m")],-0.04,"*")
text(z[which(data[,2]=="ZP562s1")],-0.04,"*")
text(z[which(data[,2]=="ZP560s1")],-0.04,"*")
#text(z[which(data[,2]=="ZP636s1")],-0.04,"*")
#text(z[which(data[,2]=="ZP565s1")],-0.04,"*")
#text(z[which(data[,2]=="ZP561s1")],-0.04,"*")
#text(z[which(data[,2]=="BC187m")],-0.04,"*")
#text(z[which(data[,2]=="L1528m")],-0.04,"*")
#text(z[which(data[,2]=="DBVPG6765m")],-0.04,"*")
#text(z[which(data[,2]=="DBVPG1788m")],-0.04,"*")
#text(z[which(data[,2]=="DBVPG1373m")],-0.04,"*")
#text(z[which(data[,2]=="YJM981m")],-0.04,"*")
#text(z[which(data[,2]=="YJM978m")],-0.04,"*")
#text(z[which(data[,2]=="YJM975m")],-0.04,"*")
#text(z[which(data[,2]=="RM11")],-0.04,"*")

# Figure C
data<-read.table("allCensFrom_gfa2.infile.k2.structure.data")
attach(data)
x<-rbind(data[,6],data[,7])
par(mar=c(7, 7, 4, 5) + 0.1,cex=0.6)
#barplot(x,beside=FALSE,col=c("#005a32","#5aae61"),names.arg=data[,2],las=2, border=NA,space=0,main='c) North American oaks, K=2',xaxs = "i")
barplot(x,beside=FALSE,col=c("#005a32","#5aae61"),names.arg=data[,2],las=2, border=NA,space=0,main='c) North American oak strains',xaxs = "i")


# show population names and numbers
par(cex=1)
text(5,0.5,"6")
text(17.5,0.5,"7")

# show ref strains
text(z[which(data[,2]=="YPS128m")],-0.04,"*")
text(z[which(data[,2]=="YPS606m")],-0.04,"*")
text(z[which(data[,2]=="SDO1s1")],-0.04,"*")
text(z[which(data[,2]=="SDO3s1")],-0.04,"*")
#text(z[which(data[,2]=="SDO2s1")],-0.04,"*")
#text(z[which(data[,2]=="SDO4s1")],-0.04,"*")
#text(z[which(data[,2]=="SDO6s1")],-0.04,"*")
#text(z[which(data[,2]=="SDO7s1")],-0.04,"*")
#text(z[which(data[,2]=="SDO8s1")],-0.04,"*")
#text(z[which(data[,2]=="SDO9s1")],-0.04,"*")
#text(z[which(data[,2]=="SM12s1")],-0.04,"*")
#text(z[which(data[,2]=="SM17s1")],-0.04,"*")
#text(z[which(data[,2]=="SM66s1")],-0.04,"*")
#text(z[which(data[,2]=="SM69s1")],-0.04,"*")

