pdf('figRaxml10ktree.pdf')
#png('figRaxml10ktree.png',width=1000,height=1000,res=200,pointsize=11)


library(ape)

rtree<-read.tree("RAxML_bipartitions.noadmixstrainsNov10k")		# read in tree

#plot.phylo(rtree,cex=0.5)			# figure out the name of the MAL node
##nodelabels()
#tiplabels()

par(xpd=T)							# allow text in the margins

									# root using all sequences except MAL
rrtree<-root(rtree,c(rtree$tip[1:24],rtree$tip[28:length(rtree$tip)]))

									# rotate the tree so root is at the base									
									# rotate so arrangement is ready for comparison to Figure 1
x<-c(rrtree$tip.label[25:27],rrtree$tip.label[28:30],rrtree$tip.label[23:24],rrtree$tip.label[31:57],rrtree$tip.label[1:22])
t<-rotateConstr(rrtree,x)
#plot.phylo(t,cex=0.5)	
#nodelabels()


									# color the clades by population
									# function from http://blog.phytools.org/2012/01/function-to-get-descendant-node-numbers.html
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}


euocolor<-"#a6dba0"					# set population colors
winecolor<-"#762a83"
paocolor<-"#005a32"
ncocolor<-"#5aae61"
malcolor<-"#377eb8"
sakecolor<-"#e41a1c"
wacolor<-"orange"
fruitcolor<-"#dd1c77"			# dark pink
unknowncolor<-"#bdbdbd" 		# grey
humancolor<-"#fa9fb5"			# pink
fermentcolor<-"black"

pew<-2								# population edge width
winenode<-99
euonode<-91
paonode<-73
nconode<-60
malnode<-85
sakenode<-87
wanode<-84

ec<-rep("black",length(t$edge[,2]))	# default edge color = black
ew<-rep(1,length(t$edge[,2]))		# default edge width = 1
tc<-rep("black",length(t$tip))		# default tip color = black (fermentations)
tc2<-rep("black",length(t$tip))		# default tip color = black (fermentations)
tc[grep("SM",t$tip)]<-ncocolor
tc[grep("SD",t$tip)]<-ncocolor
tc[grep("YPS",t$tip)]<-paocolor
tc[grep("ZP",t$tip)]<-euocolor
tc[grep("PYR",t$tip)]<-euocolor
tc[grep("HUN",t$tip)]<-euocolor
tc[grep("UWOPS0",t$tip)]<-fruitcolor
tc[grep("ARN",t$tip)]<-fruitcolor
tc[grep("AN",t$tip)]<-fruitcolor
tc[grep("DBVPG1106",t$tip)]<-fruitcolor
tc[grep("Y55",t$tip)]<-fruitcolor

tc[grep("DBVPG1788m",t$tip)]<-unknowncolor
tc[grep("DBVPG6765m",t$tip)]<-unknowncolor
tc[grep("DBVPG1373m",t$tip)]<-unknowncolor
tc[grep("YJM978m",t$tip)]<-humancolor
tc[grep("YJM975m",t$tip)]<-humancolor
tc[grep("YJM981m",t$tip)]<-humancolor

#tc[38:55]<-winecolor				# wine population
d<-getDescendants(t,node=winenode)		
for (x in d) { ec[t$edge[,2]==x]<-winecolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-winecolor }
ec[t$edge[,2]==winenode]<-winecolor		# color the edge leading to this clade
ew[t$edge[,2]==winenode]<-pew

#tc[29:37]<-euocolor					# euo population
d<-getDescendants(t,node=euonode)		
for (x in d) { ec[t$edge[,2]==x]<-euocolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-euocolor }
ec[t$edge[,2]==euonode]<-euocolor
ew[t$edge[,2]==euonode]<-pew

#tc[14:23]<-paocolor					# pao population
d<-getDescendants(t,node=paonode)	
for (x in d) { ec[t$edge[,2]==x]<-paocolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-paocolor }
ec[t$edge[,2]==paonode]<-paocolor
ew[t$edge[,2]==paonode]<-pew

#tc[c(1:13,77:80)]<-ncocolor			# nco population
d<-getDescendants(t,node=nconode)	
for (x in d) { ec[t$edge[,2]==x]<-ncocolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-ncocolor }
ec[t$edge[,2]==nconode]<-ncocolor
ew[t$edge[,2]==nconode]<-pew

									# mal population
d<-getDescendants(t,node=malnode)	
for (x in d) { ec[t$edge[,2]==x]<-malcolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-malcolor }
ec[t$edge[,2]==malnode]<-malcolor
ew[t$edge[,2]==malnode]<-pew

									# sake population
d<-getDescendants(t,node=sakenode)	
for (x in d) { ec[t$edge[,2]==x]<-sakecolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-sakecolor }
ec[t$edge[,2]==sakenode]<-sakecolor
ew[t$edge[,2]==sakenode]<-pew

									# wa population
d<-getDescendants(t,node=wanode)	
for (x in d) { ec[t$edge[,2]==x]<-wacolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-wacolor }
ec[t$edge[,2]==wanode]<-wacolor
ew[t$edge[,2]==wanode]<-pew

tc2[t$tip=="DBVPG1853m"]<-winecolor

#plot(t,cex=0.5,edge.color=ec,tip.color=tc,edge.width=ew,font=2)		# colour tips by source
plot.phylo(t,cex=0.5,edge.color=ec,tip.color=tc2,edge.width=ew,font=2)		# colour tips by clade

bs<-as.numeric(t$node)				# show bootstraps >= 70%
bs[bs<70]<-NA

# MANUALLY FIX THE BOOTSTRAPS
bs[24]<-NA
bs[14]<-96
bs[2]<-97


nodelabels(bs,frame="none",cex=0.6,bg="white",adj = c(1.1,-0.2))
add.scale.bar(0,-3,cex=0.8)

#text(0.012,55, "North Carolina Oak", col=ncocolor)
#text(0.0088,42, "Pennsylvania Oak", col=paocolor)
#text(0.0155,35, "Sake", col=sakecolor)
#text(0.0124,32.5, "West Africa", col=wacolor)
#text(0.0155,22, "Wine", col=winecolor)
#text(0.0155,9, "European Oak", col=euocolor)
#text(0.0153,2, "Malaysia", col=malcolor)


tc<-rep("black",length(t$tip))		# default tip color = black (fermentations)
tc[grep("SM",t$tip)]<-ncocolor
tc[grep("SD",t$tip)]<-ncocolor
tc[grep("YPS",t$tip)]<-paocolor
tc[grep("ZP",t$tip)]<-euocolor
tc[grep("PYR",t$tip)]<-euocolor
tc[grep("HUN",t$tip)]<-euocolor
tc[grep("UWOPS0",t$tip)]<-fruitcolor
tc[grep("ARN",t$tip)]<-fruitcolor
tc[grep("AN",t$tip)]<-fruitcolor
tc[grep("DBVPG1106",t$tip)]<-fruitcolor
tc[grep("Y55",t$tip)]<-fruitcolor

tc[grep("DBVPG1788m",t$tip)]<-unknowncolor
tc[grep("DBVPG6765m",t$tip)]<-unknowncolor
tc[grep("DBVPG1373m",t$tip)]<-unknowncolor
tc[grep("YJM978m",t$tip)]<-humancolor
tc[grep("YJM975m",t$tip)]<-humancolor
tc[grep("YJM981m",t$tip)]<-humancolor


tiplabels("",grep(ncocolor,tc),pch=20,frame="none",col=ncocolor,adj=0.5015)
tiplabels("",grep(paocolor,tc),pch=20,frame="none",col=paocolor,0.5015)
tiplabels("",grep(euocolor,tc),pch=20,frame="none",col=euocolor,0.5014)
tiplabels("",grep(fruitcolor,tc),pch=20,frame="none",col=fruitcolor,0.502)
tiplabels("",grep(unknowncolor,tc),pch=20,frame="none",col=unknowncolor,0.5021)
tiplabels("",grep(humancolor,tc),pch=20,frame="none",col=humancolor,0.5015)
tiplabels("",grep("DBVPG6044",t$tip),pch=20,frame="none",col=fermentcolor,0.5022)
tiplabels("",grep("DBVPG1853",t$tip),pch=20,frame="none",col=fermentcolor,0.5022)
tiplabels("",grep("NCYC110",t$tip),pch=20,frame="none",col=fermentcolor,0.5019)
tiplabels("",grep("L1",t$tip),pch=20,frame="none",col=fermentcolor,0.5013)
tiplabels("",grep("RM11",t$tip),pch=20,frame="none",col=fermentcolor,0.5011)
tiplabels("",grep("K11",t$tip),pch=20,frame="none",col=fermentcolor,0.5011)
tiplabels("",grep("Y12",t$tip),pch=20,frame="none",col=fermentcolor,0.5011)
tiplabels("",grep("Y9",t$tip),pch=20,frame="none",col=fermentcolor,0.5010)

#tiplabels("",c(grep(ncocolor,tc),grep("YPS13",t$tip),grep("YPS3",t$tip),grep("YPS4",t$tip),grep("YPS61",t$tip),grep("YPS600",t$tip),grep("YPS602",t$tip),grep("YPS604",t$tip),grep("YPS608",t$tip),grep("YPS610",t$tip),grep(euocolor,tc)),pch=20,frame="none",col=oakcolor,0.5011)



#tiplabels("",grep("AN",t$tip),pch=20,frame="none",col=fruitcolor,0.5012)
#tiplabels("",grep("ARN",t$tip),pch=20,frame="none",col=fruitcolor,0.5015)


lineagenames<-c("Malaysia","Sake","West Africa","Wine","North Carolina Oak","Pennsylvania Oak","European Oak")
lineagecolors<-c(malcolor,sakecolor,wacolor,winecolor,ncocolor,paocolor,euocolor)
lineagelty<-rep(1,7)
lineagelwd<-rep(2,7)
#legend(0.0105,56,lineagenames,col=lineagecolors,lty=lineagelty,lwd=lineagelwd,bty="n",cex=0.8,title="Populations",title.adj=0)
legend(0.015,58,lineagenames,col=lineagecolors,lty=lineagelty,lwd=lineagelwd,bty="n",cex=0.8,title="Populations",title.adj=0,text.col=lineagecolors,title.col="black")

sourcenames<-c("Fermentation","Fruit or flower","Human Infection","Soil or unknown","North Carolina Oak","Pennsylvania Oak","European Oak")
sourcecolors<-c("black",fruitcolor,humancolor,unknowncolor,ncocolor,paocolor,euocolor)
#legend(0.016,56,sourcenames,text.col=sourcecolors,bty="n",cex=0.8,title="Strain source")
sourcepch<-rep(20,7)
legend(0.015,4,sourcenames,pch=sourcepch,col=sourcecolors,bty="n",cex=0.8,title="Strain source",title.adj=0)


dev.off()

