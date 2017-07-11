pdf('figNJtree.pdf')

library(ape)

njtree<-read.tree("alcat.tree")		# read in tree

par(xpd=T)							# allow text in the margins


#plot.phylo(njtree,cex=0.5)			# figure out the name of the MAL node
#nodelabels(as.character((1:length(njtree$node))+length(njtree$tip)))

rnjtree<-root(njtree,node=109)		# root using the MAL node

#plot.phylo(rnjtree,cex=0.5)	
									# rotate the tree so root is at the base
x<-c(rev(rnjtree$tip.label)[5:length(rnjtree$tip.label)],rev(rnjtree$tip.label)[1:4]) # figure out a nice looking tip order
t<-rotateConstr(rnjtree,x)

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


europecolor<-"#8c510a"
usacolor<-"#4daf4a"
#euocolor<-"#a6dba0"					# set population colors
#winecolor<-"#762a83"
#paocolor<-"#005a32"
#ncocolor<-"#5aae61"
malcolor<-"#377eb8"
sakecolor<-"#e41a1c"
wacolor<-"orange"
pew<-2								# population edge width
#oakcolor<-"#543005"					# colors to highlight main samples with coloured points
#fruitcolor<-"#f46d43"

fruitcolor<-"#dd1c77"			# dark pink
unknowncolor<-"#bdbdbd" 		# grey
humancolor<-"#fa9fb5"			# pink

usanode<-97
europenode<-118

ec<-rep("black",length(t$edge[,2]))	# default edge color = black
ew<-rep(1,length(t$edge[,2]))		# default edge width = 1
tc<-rep("black",length(t$tip))		# default tip color = black (fermentations)
tc2<-rep("black",length(t$tip))		# COLOURING TIPS BY CLADE: default tip color = black (fermentations)

tc[grep("SM",t$tip)]<-ncocolor
tc[grep("SD",t$tip)]<-ncocolor
tc[grep("YPS",t$tip)]<-paocolor
tc[grep("ZP",t$tip)]<-euocolor
tc[grep("PYR",t$tip)]<-euocolor
tc[grep("HUN",t$tip)]<-euocolor
tc[grep("UWOPS",t$tip)]<-fruitcolor
tc[grep("ARN",t$tip)]<-fruitcolor
tc[grep("AN",t$tip)]<-fruitcolor
tc[grep("DBVPG1106",t$tip)]<-fruitcolor
tc[grep("Y55",t$tip)]<-fruitcolor

tc[grep("DBVPG1788m",t$tip)]<-unknowncolor
tc[grep("DBVPG6765m",t$tip)]<-unknowncolor
tc[grep("DBVPG1373m",t$tip)]<-unknowncolor
tc[grep("SK1",t$tip)]<-unknowncolor
tc[grep("YJM978m",t$tip)]<-humancolor
tc[grep("YJM975m",t$tip)]<-humancolor
tc[grep("YJM981m",t$tip)]<-humancolor
tc[grep("YJM981m",t$tip)]<-humancolor
tc[grep("273614Nm",t$tip)]<-humancolor
tc[grep("322134Sm",t$tip)]<-humancolor
tc[grep("378604Xm",t$tip)]<-humancolor
tc[grep("YJM789",t$tip)]<-humancolor
tc[grep("S288c",t$tip)]<-fruitcolor


									# europe population
d<-getDescendants(t,node=europenode)		
for (x in d) { ec[t$edge[,2]==x]<-europecolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-europecolor }
ec[t$edge[,2]==europenode]<-europecolor		# color the edge leading to this clade
ew[t$edge[,2]==europenode]<-pew


									# usa population
d<-getDescendants(t,node=usanode)	
for (x in d) { ec[t$edge[,2]==x]<-usacolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-usacolor }
ec[t$edge[,2]==usanode]<-usacolor
ew[t$edge[,2]==usanode]<-pew


									# mal population
d<-getDescendants(t,node=156)	
for (x in d) { ec[t$edge[,2]==x]<-malcolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-malcolor }
ec[t$edge[,2]==156]<-malcolor
ew[t$edge[,2]==156]<-pew

									# sake population
d<-getDescendants(t,node=151)	
for (x in d) { ec[t$edge[,2]==x]<-sakecolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-sakecolor }
ec[t$edge[,2]==151]<-sakecolor
ew[t$edge[,2]==151]<-pew

									# wa population
d<-getDescendants(t,node=155)	
for (x in d) { ec[t$edge[,2]==x]<-wacolor; ew[t$edge[,2]==x]<-pew; tc2[x]<-wacolor }
ec[t$edge[,2]==155]<-wacolor
ew[t$edge[,2]==155]<-pew

#plot(t,cex=0.5,edge.color=ec,tip.color=tc,edge.width=ew,font=2)		# colour tips by source
plot(t,cex=0.5,edge.color=ec,tip.color=tc2,edge.width=ew,font=2)		# colour tips by clade

bs<-as.numeric(t$node)				# show bootstraps >= 70%
bs[bs<70]<-NA
nodelabels(bs,frame="none",cex=0.5,bg="white",adj = c(1.1,-0.2))
add.scale.bar(0,-1)

# Figure legends

lineagenames<-c("Malaysia","Sake","West Africa","Europe","USA")
lineagecolors<-c(malcolor,sakecolor,wacolor,europecolor,usacolor)
lineagelty<-rep(1,5)
lineagelwd<-rep(2,5)
legend(0.01,80,lineagenames,col=lineagecolors,lty=lineagelty,lwd=lineagelwd,bty="n",cex=0.8,title="Populations",title.adj=0,text.col=lineagecolors,title.col="black")

#sourcenames<-c("Fermentation","Fruit or flower","Human Infection","Soil or unknown","North Carolina Oak","Pennsylvania Oak","European Oak")
#sourcecolors<-c("black",fruitcolor,humancolor,unknowncolor,ncocolor,paocolor,euocolor)
#legend(0.011,80,sourcenames,text.col=sourcecolors,bty="n",cex=0.8,title="Strain source",title.adj=0)


#tiplabels("",c(grep(ncocolor,tc),grep("YPS13",t$tip),grep("YPS3",t$tip),grep("YPS4",t$tip),grep("YPS61",t$tip),grep("YPS600",t$tip),grep("YPS602",t$tip),grep("YPS604",t$tip),grep("YPS608",t$tip),grep("YPS610",t$tip),grep(euocolor,tc)),pch=20,frame="none",col=oakcolor,0.5011)

#tiplabels("",grep("AN",t$tip),pch=20,frame="none",col=fruitcolor,0.5012)
#tiplabels("",grep("ARN",t$tip),pch=20,frame="none",col=fruitcolor,0.5015)

