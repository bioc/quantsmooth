# Support functions for SnpSetIllumina
numericCHR<- function(CHR) {
  # Set autosomal chromosomes to their number
  # X - 98
  # Y - 99
  # XY - 100
  CHR<-as.character(CHR)
  CHR[CHR=="X"]<-"98"
  CHR[CHR=="Y"]<-"99"
  CHR[CHR=="XY"]<-"100"
  as.numeric(CHR)
}
#
characterCHR<- function(CHR) {
  CHR<-as.character(CHR)
  CHR[CHR=="98"]<-"X"
  CHR[CHR=="99"]<-"Y"
  CHR[CHR=="100"]<-"XY"
  CHR
}
#
scaleto <-function(x,fromlimits=c(0,50),tolimits=c(0.5,-0.5),adjust=TRUE) {
  if (adjust) {
    x[x>fromlimits[2]]<-fromlimits[2]
    x[x<fromlimits[1]]<-fromlimits[1]
  }  
  x<- x-fromlimits[1]
  x<- x/(fromlimits[2]-fromlimits[1])
  x<- x * (tolimits[2]-tolimits[1])
  x+tolimits[1]
}  
#
plotChromosome<-function(gendata,chrompos,chromosome,dataselection=NULL,ylim=NULL,normalized.to=NULL,grid=NULL,smooth.lambda=2,interval=0.5,...) {
  # uses gcsmoothing.R
  if (is.null(dataselection)) dataselection<-rep(TRUE,ncol(gendata))
  plotSmoothed(gendata[chrompos[,"CHR"]==chromosome,dataselection],chrompos[chrompos[,"CHR"]==chromosome,"MapInfo"],ylim,normalized.to,grid,smooth.lambda,interval,...)
}
#
prepareGenomePlot<-function(chrompos,cols="grey50",paintCytobands=FALSE,bleach=0,topspace=1,organism,sexChromosomes=FALSE,...) {
  # prepare plot with lines and axes indicating all chromosomes
  # sends extra arguments to plot function
  cytobandWidth<-0.075
  # hsa 22+ XY
  # mmu 19 + XY
  # rno 20 + XY
	par(mar=c(1,4,2,3)+0.1)

	if (!missing(organism)) {
	  organism<-match.arg(organism,c("hsa","mmu","rno"))
	  chrom.n<-switch(organism,
	                   hsa = 22,
                     mmu = 19,
                     rno = 20) 
  	chrs2<-factor(numericCHR(chrompos[,"CHR"]),levels=c(1:chrom.n,if(sexChromosomes)c(98,99)else NULL))
  	lens<-sapply(split(chrompos[,"MapInfo"],chrs2),function(x)max(c(0,x)))
  	cols<-rep(cols,length.out=length(lens))
  	names(cols)<-names(lens)
  	dwidth<-NULL
  	# plot 2 columns of chromosomes, first column large->small (1-12), second column small->large (22-12)
  	for (i in 1:(chrom.n %/% 2)) dwidth[i]<-lens[i]+lens[chrom.n+1-i]
    # make sure vector length equals nr of rows in plot
  	if (chrom.n %% 2 ==1) dwidth<-c(dwidth,lens[chrom.n %/% 2 +1])
  	if (sexChromosomes) dwidth<-c(dwidth,lens["98"]+lens["99"])
  	maxdwidth<-max(dwidth)*1.05
  	leftrow<-c(if(sexChromosomes)"98" else NULL,((chrom.n + 1) %/% 2):1)
  	rightrow<-c(if(sexChromosomes)"99" else NULL, if (chrom.n %% 2 ==1) "" else NULL,((chrom.n + 1) %/% 2 +1):chrom.n)
  	plot(c(0,maxdwidth),c(0.5 ,0.5+length(dwidth)+topspace),type="n",ylab="Chromosome",xlab="",axes = FALSE, las = 2,...)
  	axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
  	axis(4, c(1:length(dwidth)), characterCHR(rightrow), las = 2)
  	if (paintCytobands && organism=="hsa") {
    	for (i in 1:length(dwidth)) {
    	  if (lens[leftrow[i]]>0) paintCytobands(i,c(0,i+cytobandWidth/2),"bases",width=cytobandWidth,length.out=lens[leftrow[i]],legend=FALSE,bleach=bleach)
    	  if(rightrow[i]!="" && lens[rightrow[i]]>0) paintCytobands(i,c(maxdwidth-lens[rightrow[i]],i+cytobandWidth/2),"bases",width=cytobandWidth,length.out=lens[rightrow[i]],legend=FALSE,bleach=bleach)
  	  }
  	} else {
    	for (i in 1:length(dwidth)) {
    	  lines(c(0,lens[leftrow[i]]),c(i,i),col=cols[leftrow[i]],lwd=2)
    	  if(rightrow[i]!="") lines(c(maxdwidth-lens[rightrow[i]],maxdwidth),c(i,i),col=cols[rightrow[i]],lwd=2)
    	}
    }
    # for each locus determine postion on plot , this can be used later to fill with data
  	dchrompos<-matrix(0,nrow=nrow(chrompos),ncol=2,dimnames=list(rownames(chrompos),c("CHR","MapInfo")))
 		for (i in 1:length(rightrow)) if (rightrow[i]!="") {
 		  probes<-numericCHR(chrompos[,"CHR"])==rightrow[i]
      dchrompos[probes,2]<-chrompos[probes,"MapInfo"]+maxdwidth-lens[rightrow[i]]
      dchrompos[probes,1]<- i
    }
  	for (i in 1:length(leftrow)) {
 		  probes<-numericCHR(chrompos[,"CHR"])==leftrow[i]
			dchrompos[probes,2]<-chrompos[probes,"MapInfo"]
      dchrompos[probes,1]<- i
    }
	}
  else {
  	chrs2<-factor(numericCHR(chrompos[,"CHR"]))
  	lens<-sapply(split(chrompos[,"MapInfo"],chrs2),max)
  	m<-length(lens)
  	cols<-rep(cols,length.out=m)
    maxdwidth<-max(lens)
  	plot(c(0,maxdwidth),c(0.5,m+0.5+topspace),type="n",ylab="Chromosome",xlab="",axes = FALSE, las = 2,...)
  	axis(2, c(m:1), characterCHR(names(lens)), las = 2)
    for (i in 1:m)  lines(c(0,lens[i]),c(m+1-i,m+1-i),col=cols[as.numeric(names(lens))],lwd=2)
    dchrompos<-chrompos
    dchrompos[,1]<-m+1-as.numeric(chrs2)
  }
	dchrompos
}
#  data taken from lodplot package
#  original data available at: ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/mapview/BUILD.35.1/ideogram.gz.
semicircle <- function(base.x, base.y, base.length, height=base.length, side=1, orientation=NULL,plottype="poly",...) {
  # based on lodplot package
  # - col is now propagated through ..., other plotting parameters can now also be given
  # - different types poly/line 
  radius<-base.length/2
  x<-radius*seq(-1,1,length=40)
  y<-height/radius*sqrt(radius^2-x^2)
  if (is.null(orientation)) {
    co<-as.integer(cos(pi*(3-side)/2))
    so<-as.integer(sin(pi*(3-side)/2))
  }else{
    co<-cos(orientation)
    so<-sin(orientation)
  }
  tx<-co*x - so*y
  ty<-so*x + co*y
  if (is.null(orientation)) {
    if (side==1 || side==3) {
      base.x<-base.x+radius
    }else if (side==2 || side==4) {
      base.y<-base.y+radius
    }
  }
  x<-base.x+tx
  y<-base.y+ty
  switch(plottype,
    poly=polygon(x,y,...),
    line=lines(x,y,...)
  )
}
#
paintCytobands<-function(chrom, pos=c(0,0), units=c("cM","bases","ISCN"), width=0.4, length.out, bands="major", orientation=c("h","v"), legend = TRUE, cex.leg=0.7, bleach = 0) {
  # Based on paint.chromosome from lodplot package
  # added:
  #  -bleach
  #  -length.out
  #  -using all of cM,bases,ISCN
  #  -using hatches for stalk, acen
  #  -legend + cex.leg
  #  -orientation
  #  extracted semicircle for general use
  bleacher<-function(x) { (x * (1-bleach)) + bleach}
  require(lodplot)
  data(chrom.bands)
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  units<-match.arg(units)
  orientation<-match.arg(orientation)
  # original function only required ypos
  if (length(pos)==1) pos<-c(0,pos)
  chromdata<-subset(chrom.bands, chrom.bands$chr==chrom)
  lc<-nchar(chromdata$band)
  sel<-!(substr(chromdata$band,lc,lc) %in% letters)
  if (bands!="major") sel<-!sel
  chromdata<-chromdata[sel,]
  rm(lc,sel)
  bandpos<-switch(units,
         cM =chromdata[,c("cM.top","cM.bot")],
         bases = chromdata[,c("bases.top","bases.bot")],
         ISCN =  chromdata[,c("ISCN.top","ISCN.bot")])

  type.b<-match(chromdata$stain,c("acen","gneg", "gpos", "gvar", "stalk"))
  bandcol<-gray(bleacher(c(0.5,1,0.2,0.6,0.75)))[type.b]
  banddens<-c(30,-1,-1,-1,10)[type.b]
  bandbord<-gray(bleacher(c(0,0,0,0,1)))[type.b]
  if (!missing(length.out)) {
    bandpos<-(bandpos/max(bandpos))*length.out
  }                                          
  n<-nrow(chromdata)
  centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
  idx<-c(2:(centromere-1), (centromere+2):(n-1))
  if (orientation=="h") {
    rect(pos[1]+bandpos[idx,1],pos[2],pos[1]+bandpos[idx,2],pos[2]-width, col=bandcol[idx], density=banddens[idx], border=bandbord[idx])
    semicircle(pos[1]+bandpos[1,2], pos[2]-width, width,
               bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=bandbord[1])
    semicircle(pos[1]+bandpos[n,1], pos[2]-width, width,
               bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=bandbord[n])
    semicircle(pos[1]+bandpos[centromere,1], pos[2]-width, width,
               bandpos[centromere,2]-bandpos[centromere,1],
               4, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere])
    semicircle(pos[1]+bandpos[centromere+1,2], pos[2]-width, width,
               bandpos[centromere+1,2]-bandpos[centromere+1,1],
               2, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1])

    centromere.size=0.6*0.5*width/yinch(1)
    symbols(pos[1]+bandpos[centromere,2], pos[2]-0.5*width,circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white")
    if (legend) text(pos[1]+(bandpos[,1]+bandpos[,2])/2,pos[2]+0.5*width,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=90,cex=cex.leg)
  } else {
    rect(pos[1],pos[2]-bandpos[idx,1],pos[1]-width,pos[2]-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=bandbord[idx])
    semicircle(pos[1]-width, pos[2]-bandpos[1,2], width,
               bandpos[1,2]-bandpos[1,1], 3, col=bandcol[1], density=banddens[1], border=bandbord[1])
    semicircle(pos[1]-width, pos[2]-bandpos[n,1], width,
               bandpos[n,2]-bandpos[n,1], 1, col=bandcol[n], density=banddens[n], border=bandbord[n])
    semicircle(pos[1]-width, pos[2]-bandpos[centromere,1], width,
               bandpos[centromere,2]-bandpos[centromere,1],
               1, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere])
    semicircle(pos[1]-width, pos[2]-bandpos[centromere+1,2], width,
               bandpos[centromere+1,2]-bandpos[centromere+1,1],
               3, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1])
    centromere.size=0.6*0.5*width/xinch(1)
    symbols(pos[1]-0.5*width, pos[2]-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white")
    if (legend) text(pos[1]+0.5*width,pos[2]-(bandpos[,1]+bandpos[,2])/2,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=0,cex=cex.leg)
  }
}

lengthChromosome<-function(chrom, units=c("cM","bases","ISCN")) {
  require(lodplot)
  data(chrom.bands)
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  units<-match.arg(units)
  chromdata<-subset(chrom.bands, chrom.bands$chr==chrom)
  switch(units,cM =chromdata[nrow(chromdata),"cM.bot"],
               bases = chromdata[nrow(chromdata),"bases.bot"],
               ISCN =  chromdata[nrow(chromdata),"ISCN.bot"])
}

drawSimpleChrom<-function(x,y,len=3,width=1,fill,col,orientation=c("h","v"),centromere.size=0.6) {
  # put a simple drawing of a chromosome p:q = 1:2
  # events can be indictaed by fill and col fill=c("a","p","q","p1","p2","p3","q1","q2","q3")
  bandpos<-cbind(c(0,1,2,3,4,7),c(1,2,3,4,7,8))*len/8
  n<-nrow(bandpos)
  centromere<-3
  idx<-c(2:(centromere-1), (centromere+2):(n-1))
  bandcol=rep("white",6)
  if (!missing(fill)) if(length(fill)>0) for (i in 1:length(fill)) {
    if (fill[i]=="a") bandcol[1:6]<-col[i]
    else if (fill[i]=="p") bandcol[1:3]<-col[i]
    else if (fill[i]=="q") bandcol[4:6]<-col[i]
    else if (fill[i]=="p1") bandcol[3]<-col[i]
    else if (fill[i]=="p2") bandcol[2]<-col[i]
    else if (fill[i]=="p3") bandcol[1]<-col[i]
    else if (fill[i]=="q1") bandcol[4]<-col[i]
    else if (fill[i]=="q2") bandcol[5]<-col[i]
    else if (fill[i]=="q3") bandcol[6]<-col[i]
  }  
  banddens=rep(-1,6)
  if (orientation[1]=="h") {
    # draw the inside filling
    rect(x+bandpos[idx,1],y+0.5*width,x+bandpos[idx,2],y-0.5*width, col=bandcol[idx], density=banddens[idx], border=NA)
    semicircle(x+bandpos[1,2], y-0.5*width, width,
               bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=NA)
    semicircle(x+bandpos[n,1], y-0.5*width, width,
               bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=NA)
    semicircle(x+bandpos[centromere,1], y-0.5*width, width,
               bandpos[centromere,2]-bandpos[centromere,1],
               4, col=bandcol[centromere], density=banddens[centromere], border=NA)
    semicircle(x+bandpos[centromere+1,2], y-0.5*width, width,
               bandpos[centromere+1,2]-bandpos[centromere+1,1],
               2, col=bandcol[centromere+1], density=banddens[centromere+1], border=NA)
    # draw the circumference
    for (i in idx) {
      lines(x+bandpos[i,1:2],rep(y+0.5*width,2),col=1)
      lines(x+bandpos[i,1:2],rep(y-0.5*width,2),col=1)
    }
    semicircle(x+bandpos[1,2], y-0.5*width, width, bandpos[1,2]-bandpos[1,1], 2, col=1, plottype="line")
    semicircle(x+bandpos[n,1], y-0.5*width, width, bandpos[n,2]-bandpos[n,1], 4, col=1, plottype="line")
    semicircle(x+bandpos[centromere,1], y-0.5*width, width, bandpos[centromere,2]-bandpos[centromere,1], 4, col=1, plottype="line")
    semicircle(x+bandpos[centromere+1,2], y-0.5*width, width, bandpos[centromere+1,2]-bandpos[centromere+1,1], 2, col=1, plottype="line")
    # draw the centromere
    centromere.size=centromere.size*0.5*width/yinch(1)
    symbols(x+bandpos[centromere,2], y,circles=1,inches=centromere.size, add=TRUE,fg="black",bg="white")
  } else {
    # draw the inside filling
    rect(x+0.5*width,y-bandpos[idx,1],x-0.5*width,y-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=NA)
    semicircle(x-0.5*width, y-bandpos[1,2], width, bandpos[1,2]-bandpos[1,1],
               3, col=bandcol[1], density=banddens[1], border=NA)
    semicircle(x-0.5*width, y-bandpos[n,1], width, bandpos[n,2]-bandpos[n,1], 
               1, col=bandcol[n], density=banddens[n], border=NA)
    semicircle(x-0.5*width, y-bandpos[centromere,1], width, bandpos[centromere,2]-bandpos[centromere,1],
               1, col=bandcol[centromere], density=banddens[centromere], border=NA)
    semicircle(x-0.5*width, y-bandpos[centromere+1,2], width, bandpos[centromere+1,2]-bandpos[centromere+1,1],
               3, col=bandcol[centromere+1], density=banddens[centromere+1], border=NA)
    # draw the circumference
    for (i in idx) {
      lines(rep(x+0.5*width,2),y-bandpos[i,1:2],col=1)
      lines(rep(x-0.5*width,2),y-bandpos[i,1:2],col=1)
    }
    semicircle(x-0.5*width, y-bandpos[1,2], width, bandpos[1,2]-bandpos[1,1], 3, col=1, plottype="line")
    semicircle(x-0.5*width, y-bandpos[n,1], width, bandpos[n,2]-bandpos[n,1], 1, col=1, plottype="line")
    semicircle(x-0.5*width, y-bandpos[centromere,1], width, bandpos[centromere,2]-bandpos[centromere,1], 1, col=1, plottype="line")
    semicircle(x-0.5*width, y-bandpos[centromere+1,2], width, bandpos[centromere+1,2]-bandpos[centromere+1,1], 3, col=1, plottype="line")
    # draw the centromere
    centromere.size=centromere.size*0.5*width/xinch(1)
    symbols(x, y-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg="black",bg="white")
  }
}



