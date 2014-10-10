# Support functions for SnpSetIllumina
numericCHR<- function(CHR, prefix="chr") {
  # Set autosomal chromosomes to their number
  # X - 98
  # Y - 99
  # XY - 100
  # MT, M - 101
  # Variable regions - 102
  CHR<-sub(paste0("^",prefix),"",as.character(CHR))
  CHR[grep("_",CHR)]<-"102"
  CHR[CHR=="X"]<-"98"
  CHR[CHR=="Y"]<-"99"
  CHR[CHR=="XY"]<-"100"
  CHR[CHR %in% c("M","MT")]<-"101"
  as.numeric(CHR)
}
#
characterCHR<- function(CHR,prefix="") {
  CHR<-as.character(CHR)
  CHR[CHR=="98"]<-"X"
  CHR[CHR=="99"]<-"Y"
  CHR[CHR=="100"]<-"XY"
  CHR[CHR=="101"]<-"MT"
  CHR[CHR=="102"]<-"-" 
  paste(prefix,CHR,sep="")
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
.getChrombands<-function(units) {
  if (units %in% c("cM","bases","ISCN")) {
    chrom.bands<-NULL;rm(chrom.bands) # trick to satisfy R check
    data(chrom.bands,package="quantsmooth",envir=environment())
    bandpos<-switch(units,
                    cM =chrom.bands[,c("cM.top","cM.bot")],
                    bases = chrom.bands[,c("bases.top","bases.bot")],
                    ISCN =  chrom.bands[,c("ISCN.top","ISCN.bot")])
    data.frame(chr=chrom.bands$chr,segstart=bandpos[,1],segend=bandpos[,2],stain=chrom.bands$stain,band=chrom.bands$band,arm=chrom.bands$arm, stringsAsFactors=FALSE)
  } else {
    data(list=paste0("chrom.bands.",units),package="quantsmooth",envir=environment())
    .convertUCSCcytoband(get(paste0("chrom.bands.",units)))
  }
}
#
.convertUCSCcytoband<- function(ucscdata) {
  data.frame(chr=characterCHR(numericCHR(ucscdata[,1])),segstart=ucscdata[,2],segend=ucscdata[,3],stain=ucscdata[,5],band=substring(ucscdata[,4],2),arm=substring(ucscdata[,4],1,1), stringsAsFactors=FALSE)
}
#
plotChromosome<-function(gendata,chrompos,chromosome,dataselection=NULL,ylim=NULL,normalized.to=NULL,grid=NULL,smooth.lambda=2,interval=0.5,...) {
  # uses gcsmoothing.R
  if (is.null(dataselection)) dataselection<-rep(TRUE,ncol(gendata))
  plotSmoothed(gendata[chrompos[,"CHR"]==chromosome,dataselection],chrompos[chrompos[,"CHR"]==chromosome,"MapInfo"],ylim=ylim,normalized.to=normalized.to,
               grid=grid,smooth.lambda=smooth.lambda,interval=interval,...)
}
#
prepareGenomePlot<-function(chrompos=NULL,cols="grey50",paintCytobands=FALSE,bleach=0,topspace=1,organism,sexChromosomes=FALSE,units="hg19",...) {
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
    if (is.null(chrompos)) {
      chroms=c(1:chrom.n,if(sexChromosomes)c(98,99)else NULL)
      chrnames=characterCHR(chroms,"chr")
      mapinfo=rep(0,length(chroms))
      chrompos=cbind(CHR=chroms,MapInfo=mapinfo)
      rownames(chrompos)<-chrnames
    }
  	chrs2<-factor(numericCHR(chrompos[,"CHR"]),levels=c(1:chrom.n,if(sexChromosomes)c(98,99)else NULL))
  	if (organism=="hsa")
      lens<-lengthChromosome(levels(chrs2),units=units)
  	else
    	lens<-sapply(split(chrompos[,"MapInfo"],chrs2),function(x)max(c(0,x)))
  	names(lens)<-characterCHR(names(lens))
  	cols<-rep(cols,length.out=length(lens))
  	names(cols)<-names(lens)
  	dwidth<-NULL
  	# plot 2 columns of chromosomes, first column large->small (1-12), second column small->large (22-12)
  	for (i in 1:(chrom.n %/% 2)) dwidth[i]<-lens[i]+lens[chrom.n+1-i]
    # make sure vector length equals nr of rows in plot
  	if (chrom.n %% 2 ==1) dwidth<-c(dwidth,lens[chrom.n %/% 2 +1])
  	if (sexChromosomes) dwidth<-c(dwidth,lens["X"]+lens["Y"])
  	maxdwidth<-max(dwidth)*1.05
  	leftrow<-c(if(sexChromosomes)"X" else NULL,((chrom.n + 1) %/% 2):1)
  	rightrow<-c(if(sexChromosomes)"Y" else NULL, if (chrom.n %% 2 ==1) "" else NULL,((chrom.n + 1) %/% 2 +1):chrom.n)
  	plot(c(0,maxdwidth),c(0.5 ,0.5+length(dwidth)+topspace),type="n",ylab="Chromosome",xlab="",axes = FALSE, las = 2,...)
  	axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
  	axis(4, c(1:length(dwidth)), characterCHR(rightrow), las = 2)
  	if (paintCytobands && organism=="hsa") {
    	for (i in 1:length(dwidth)) {
    	  if (lens[leftrow[i]]>0) paintCytobands(leftrow[i],c(0,i+cytobandWidth/2),units,width=cytobandWidth,length.out=lens[leftrow[i]],legend=FALSE,bleach=bleach)
    	  if (rightrow[i]!="" && lens[rightrow[i]]>0) paintCytobands(rightrow[i],c(maxdwidth-lens[rightrow[i]],i+cytobandWidth/2),"bases",width=cytobandWidth,length.out=lens[rightrow[i]],legend=FALSE,bleach=bleach)
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
 		  probes<-characterCHR(chrompos[,"CHR"])==rightrow[i]
      dchrompos[probes,2]<-chrompos[probes,"MapInfo"]+maxdwidth-lens[rightrow[i]]
      dchrompos[probes,1]<- i
    }
  	for (i in 1:length(leftrow)) {
 		  probes<-characterCHR(chrompos[,"CHR"])==leftrow[i]
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
qs.semicircle <- function(base.x, base.y, base.length, height=base.length, side=1, orientation=NULL,plottype="poly",...) {
  # based on lodplot package
  # David Duffy <David.Duffy@qimr.edu.au>
  # URL: http://www.qimr.edu.au/davidD
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

qs.grid.semicircle<-function (base.x, base.y, base.length, height = base.length, 
    side = 1, orientation = NULL, ...) 
{
    radius <- base.length/2
    x <- radius * seq(-1, 1, length = 40)
    y <- height/radius * sqrt(radius^2 - x^2)
    if (is.null(orientation)) {
        co <- as.integer(cos(pi * (3 - side)/2))
        so <- as.integer(sin(pi * (3 - side)/2))
    }
    else {
        co <- cos(orientation)
        so <- sin(orientation)
    }
    tx <- co * x - so * y
    ty <- so * x + co * y
    if (is.null(orientation)) {
        if (side == 1 || side == 3) {
            base.x <- base.x + radius
        }
        else if (side == 2 || side == 4) {
            base.y <- base.y + radius
        }
    }
    x <- base.x + tx
    y <- base.y + ty
    grid.polygon(x, y, ...)
}

#
paintCytobands<-function(chrom, pos=c(0,0), units="hg19", width=0.4, length.out, bands="major", orientation=c("h","v"), legend = TRUE, cex.leg=0.7, bleach = 0,...) {
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
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  orientation<-match.arg(orientation)
  # original function only required ypos
  if (length(pos)==1) pos<-c(0,pos)
  chromdata<-subset(chrombands, chrombands$chr==chrom)
  if (nrow(chromdata)>0){
    lc<-nchar(chromdata$band)
    sel<-!(substr(chromdata$band,lc,lc) %in% letters)
    if (bands!="major") sel<-!sel
    chromdata<-chromdata[sel,]
    rm(lc,sel)

    type.b<-match(chromdata$stain,c("acen","gneg", "gpos", "gvar", "stalk","gpos25","gpos50","gpos75","gpos100"))
    bandpos<-chromdata[,c("segstart","segend")]
    bandcol<-gray(bleacher(c(0.5,1,0.2,0.6,0.75,0.7,0.5,0.3,0.1)))[type.b]
    banddens<-c(30,-1,-1,-1,10,-1,-1,-1,-1)[type.b]
    bandbord<-gray(bleacher(c(0,0,0,0,1,0,0,0,0)))[type.b]
    if (!missing(length.out)) {
      bandpos<-(bandpos/max(bandpos))*length.out
    }
    n<-nrow(chromdata)
    centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
    idx<-c(2:(centromere-1), (centromere+2):(n-1))
    if (orientation=="h") {
      rect(pos[1]+bandpos[idx,1],pos[2],pos[1]+bandpos[idx,2],pos[2]-width, col=bandcol[idx], density=banddens[idx], border=bandbord[idx])
      qs.semicircle(pos[1]+bandpos[1,2], pos[2]-width, width,
                 bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      qs.semicircle(pos[1]+bandpos[n,1], pos[2]-width, width,
                 bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      qs.semicircle(pos[1]+bandpos[centromere,1], pos[2]-width, width,
                 bandpos[centromere,2]-bandpos[centromere,1],
                 4, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
      qs.semicircle(pos[1]+bandpos[centromere+1,2], pos[2]-width, width,
                 bandpos[centromere+1,2]-bandpos[centromere+1,1],
                 2, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)

      centromere.size=0.6*0.5*width/yinch(1)
      symbols(pos[1]+bandpos[centromere,2], pos[2]-0.5*width,circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
      if (legend) text(pos[1]+(bandpos[,1]+bandpos[,2])/2,pos[2]+0.5*width,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=90,cex=cex.leg,...)
    } else {
      rect(pos[1],pos[2]-bandpos[idx,1],pos[1]-width,pos[2]-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=bandbord[idx],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[1,2], width,
                 bandpos[1,2]-bandpos[1,1], 3, col=bandcol[1], density=banddens[1], border=bandbord[1],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[n,1], width,
                 bandpos[n,2]-bandpos[n,1], 1, col=bandcol[n], density=banddens[n], border=bandbord[n],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[centromere,1], width,
                 bandpos[centromere,2]-bandpos[centromere,1],
                 1, col=bandcol[centromere], density=banddens[centromere], border=bandbord[centromere],...)
      qs.semicircle(pos[1]-width, pos[2]-bandpos[centromere+1,2], width,
                 bandpos[centromere+1,2]-bandpos[centromere+1,1],
                 3, col=bandcol[centromere+1], density=banddens[centromere+1], border=bandbord[centromere+1],...)
      centromere.size=0.6*0.5*width/xinch(1)
      symbols(pos[1]-0.5*width, pos[2]-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg=gray(bleacher(0)),bg="white",...)
      if (legend) text(pos[1]+0.5*width,pos[2]-(bandpos[,1]+bandpos[,2])/2,paste(chromdata[,"arm"],chromdata[,"band"],sep=""),adj=c(0,0.5),srt=0,cex=cex.leg,...)
    }
  } else {
    warning(paste("Chromosome",chrom,"is not plotted because cytoband data is not available"))
  }
}

grid.chromosome<-function (chrom, side=1, units="hg19", chrom.width=0.5, length.out, 
                           bands="major", legend = c("chrom","band","none"), cex.leg=0.7, 
                           bleach = 0, ...)
{
  bleacher<-function(x) { (x * (1-bleach)) + bleach}
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  side<-max(1,min(side,4))
  legend<-match.arg(legend)
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
    #if (new)
    #    grid.newpage()
  chromdata <- subset(chrombands, chrombands$chr == chrom)
  if (nrow(chromdata)>0){
    lc <- nchar(chromdata$band)
    sel <- !(substr(chromdata$band, lc, lc) %in% letters)
    if (bands != "major")
        sel <- !sel
    chromdata <- chromdata[sel, ]
    rm(lc, sel)
    type.b<-match(chromdata$stain,c("acen","gneg", "gpos", "gvar", "stalk","gpos25","gpos50","gpos75","gpos100"))
    bandpos<-chromdata[,c("segstart","segend")]
    bandcol<-gray(bleacher(c(0.5,1,0.2,0.6,0.75,0.7,0.5,0.3,0.1)))[type.b]
    banddens<-c(30,-1,-1,-1,10,-1,-1,-1,-1)[type.b]
    bandbord<-gray(bleacher(c(0,0,0,0,1,0,0,0,0)))[type.b]
    if (!missing(length.out)) {
      bandpos<-(bandpos/max(bandpos))*length.out
    }
    n<-nrow(chromdata)
    centromere<-which(chromdata$arm[-n]!=chromdata$arm[-1])
    idx<-c(2:(centromere-1), (centromere+2):(n-1))
    if (side %in% 1:2) {
      pos.bottom<-0
      pos.top<-chrom.width
      pos.chrom<-(1+chrom.width)/2
      pos.band<-chrom.width+0.1*(1-chrom.width)
    } else {                                                        
      pos.bottom<-1-chrom.width
      pos.top<-1
      pos.chrom<-(1-chrom.width)/2
      pos.band<-0.1*(1-chrom.width)
    }
    bleachblack<- gray(bleacher(0))
    if (side %in% c(1,3)) {
      pushViewport(viewport(xscale = c(bandpos[1,1] , bandpos[n,2] ), yscale = c(0, 1), clip = "on",...))
      grid.rect(x = bandpos[idx,1], y = pos.bottom, width = bandpos[idx,2] -
          bandpos[idx,1], height = chrom.width, just = c("left",
          "bottom"), default.units = "native", gp = gpar(fill = bandcol[idx],col=bandbord[idx]))
      qs.grid.semicircle(bandpos[1,2], pos.bottom, chrom.width,
          bandpos[1,2] - bandpos[1,1], 2, default.units="native", gp=gpar(fill = bandcol[1],col=bandbord[1]))
      qs.grid.semicircle(bandpos[n,1], pos.bottom, chrom.width,
          bandpos[n,2] - bandpos[n,1], 4, default.units="native", gp=gpar(fill = bandcol[n],col=bandbord[n]))
      qs.grid.semicircle(bandpos[centromere,1], pos.bottom, chrom.width, 
          bandpos[centromere,2] - bandpos[centromere,1],
          4, default.units="native", gp=gpar(fill = bandcol[centromere],col=bandbord[centromere]))
      qs.grid.semicircle(bandpos[centromere + 1,2], pos.bottom, chrom.width, 
          bandpos[centromere + 1,2] - bandpos[centromere + 1,1], 
          2, default.units="native", gp=gpar(fill = bandcol[centromere+1],col=bandbord[centromere+1]))
      grid.circle(bandpos[centromere,2], pos.bottom+chrom.width/2, unit(chrom.width*0.3,"npc"), default.units="native", gp = gpar(col=bleachblack, fill="white", lwd=2))
      if (legend=="chrom") {
        grid.text(chrom, unit(0.5, "npc"), unit(pos.chrom,"native"), gp = gpar(cex = cex.leg))
      } else if (legend=="band") {
        grid.text(paste(chromdata[,"arm"],chromdata[,"band"],sep=""),(bandpos[,1]+bandpos[,2])/2,pos.band,default.units="native",hjust=0,vjust=0.5,rot=90,gp=gpar(cex=cex.leg))    
      }
    } else {
      pushViewport(viewport(xscale = c(0, 1), yscale = c(bandpos[n,2], bandpos[1,1] ), clip = "on",...))
      grid.rect(x = pos.bottom, y = bandpos[idx,1], width = chrom.width, height = bandpos[idx,2] -
          bandpos[idx,1], just = c("left", "bottom"), default.units = "native", gp = gpar(fill = bandcol[idx],col=bandbord[idx]))
      qs.grid.semicircle( pos.bottom, bandpos[1,2],chrom.width,
          bandpos[1,2] - bandpos[1,1],  1, default.units="native", gp=gpar(fill = bandcol[1],col=bandbord[1]))
      qs.grid.semicircle( pos.bottom, bandpos[n,1], chrom.width,
          bandpos[n,2] - bandpos[n,1], 3, default.units="native", gp=gpar(fill = bandcol[n],col=bandbord[n]))
      qs.grid.semicircle( pos.bottom, bandpos[centromere,1],  chrom.width,
          bandpos[centromere,2] - bandpos[centromere,1], 
          3, default.units="native", gp=gpar(fill = bandcol[centromere],col=bandbord[centromere]))
      qs.grid.semicircle( pos.bottom, bandpos[centromere + 1,2],      chrom.width,
          bandpos[centromere + 1,2] - bandpos[centromere + 1,1],  
          1, default.units="native", gp=gpar(fill = bandcol[centromere+1],col=bandbord[centromere+1]))
      grid.circle(pos.bottom+chrom.width/2, bandpos[centromere,2],  unit(chrom.width*0.3,"npc"), default.units="native", gp = gpar(col=bleachblack, fill="white", lwd=2))
      if (legend=="chrom") {
        grid.text(chrom, unit(pos.chrom,"native"), unit(0.5, "npc"), gp = gpar(cex = cex.leg))
      } else if (legend=="band") {
        grid.text(paste(chromdata[,"arm"],chromdata[,"band"],sep=""),pos.band,(bandpos[,1]+bandpos[,2])/2,default.units="native",hjust=0,vjust=0.5,rot=0,gp=gpar(cex=cex.leg))    
      }
    
    }
    popViewport()
  } else {
    warning(paste("Chromosome",chrom,"is not plotted because cytoband data is not available"))
  
  }
}

lengthChromosome<-function(chrom, units="hg19") {
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  chrom<-characterCHR(chrom)
  chromdata<-subset(chrombands, chrombands$chr %in% chrom)
  if (nrow(chromdata)==0) {
    warning("chromosomes not found in chromosomal data")
    res=rep(NA,length(chrom))
  }else{
    chromlengths<-aggregate(chromdata$segend,list(chromdata$chr),function(x) x[length(x)])
    res<-chromlengths[match(chrom,chromlengths[,1]),2]
  }
  names(res)<-chrom
  return(res)
}


position2Cytoband<-function(chrom,position,units="hg19",bands=c("major","minor")) {
  if (class(units)=="data.frame") {
    chrombands<-.convertUCSCcytoband(units)
  } else{
    chrombands<-.getChrombands(units)
  }
  chrom<-switch(as.character(chrom),
         "98"="X",
         "99"="Y",
         as.character(chrom))
  bands<-match.arg(bands)
  chromdata<-subset(chrombands, chrombands$chr==chrom)
  if (nrow(chromdata)==0) stop("invalid chromosome:",chrom)
  lc<-nchar(chromdata$band)
  sel<-!(substr(chromdata$band,lc,lc) %in% letters)
  if (bands=="minor") 
    sel<-!sel
  chromdata<-chromdata[sel,]
  rm(lc,sel)
  res<-NULL
  for (pos1 in position) {         
    cb<-which(pos1>=chromdata$segstart & pos1<=chromdata$segend)
    if (length(cb)>0)         
      res<-c(res,paste(chromdata[cb,c(1,6,5)],collapse=""))
    else {
      warning("Position ",pos1," not valid for chromosome ",chrom,". It should be between ",chromdata$segstart[1]," and ",chromdata$segend[nrow(chromdata)])
      res<-c(res,"-")
    }
  }
  res
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
    qs.semicircle(x+bandpos[1,2], y-0.5*width, width,
               bandpos[1,2]-bandpos[1,1], 2, col=bandcol[1], density=banddens[1], border=NA)
    qs.semicircle(x+bandpos[n,1], y-0.5*width, width,
               bandpos[n,2]-bandpos[n,1], 4, col=bandcol[n], density=banddens[n], border=NA)
    qs.semicircle(x+bandpos[centromere,1], y-0.5*width, width,
               bandpos[centromere,2]-bandpos[centromere,1],
               4, col=bandcol[centromere], density=banddens[centromere], border=NA)
    qs.semicircle(x+bandpos[centromere+1,2], y-0.5*width, width,
               bandpos[centromere+1,2]-bandpos[centromere+1,1],
               2, col=bandcol[centromere+1], density=banddens[centromere+1], border=NA)
    # draw the circumference
    for (i in idx) {
      lines(x+bandpos[i,1:2],rep(y+0.5*width,2),col=1)
      lines(x+bandpos[i,1:2],rep(y-0.5*width,2),col=1)
    }
    qs.semicircle(x+bandpos[1,2], y-0.5*width, width, bandpos[1,2]-bandpos[1,1], 2, col=1, plottype="line")
    qs.semicircle(x+bandpos[n,1], y-0.5*width, width, bandpos[n,2]-bandpos[n,1], 4, col=1, plottype="line")
    qs.semicircle(x+bandpos[centromere,1], y-0.5*width, width, bandpos[centromere,2]-bandpos[centromere,1], 4, col=1, plottype="line")
    qs.semicircle(x+bandpos[centromere+1,2], y-0.5*width, width, bandpos[centromere+1,2]-bandpos[centromere+1,1], 2, col=1, plottype="line")
    # draw the centromere
    centromere.size=centromere.size*0.5*width/yinch(1)
    symbols(x+bandpos[centromere,2], y,circles=1,inches=centromere.size, add=TRUE,fg="black",bg="white")
  } else {
    # draw the inside filling
    rect(x+0.5*width,y-bandpos[idx,1],x-0.5*width,y-bandpos[idx,2], col=bandcol[idx], density=banddens[idx], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[1,2], width, bandpos[1,2]-bandpos[1,1],
               3, col=bandcol[1], density=banddens[1], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[n,1], width, bandpos[n,2]-bandpos[n,1], 
               1, col=bandcol[n], density=banddens[n], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[centromere,1], width, bandpos[centromere,2]-bandpos[centromere,1],
               1, col=bandcol[centromere], density=banddens[centromere], border=NA)
    qs.semicircle(x-0.5*width, y-bandpos[centromere+1,2], width, bandpos[centromere+1,2]-bandpos[centromere+1,1],
               3, col=bandcol[centromere+1], density=banddens[centromere+1], border=NA)
    # draw the circumference
    for (i in idx) {
      lines(rep(x+0.5*width,2),y-bandpos[i,1:2],col=1)
      lines(rep(x-0.5*width,2),y-bandpos[i,1:2],col=1)
    }
    qs.semicircle(x-0.5*width, y-bandpos[1,2], width, bandpos[1,2]-bandpos[1,1], 3, col=1, plottype="line")
    qs.semicircle(x-0.5*width, y-bandpos[n,1], width, bandpos[n,2]-bandpos[n,1], 1, col=1, plottype="line")
    qs.semicircle(x-0.5*width, y-bandpos[centromere,1], width, bandpos[centromere,2]-bandpos[centromere,1], 1, col=1, plottype="line")
    qs.semicircle(x-0.5*width, y-bandpos[centromere+1,2], width, bandpos[centromere+1,2]-bandpos[centromere+1,1], 3, col=1, plottype="line")
    # draw the centromere
    centromere.size=centromere.size*0.5*width/xinch(1)
    symbols(x, y-bandpos[centromere,2],circles=1,inches=centromere.size, add=TRUE,fg="black",bg="white")
  }
}



