#
# quantsmooth.r
# (c) J. Oosting 2006
.quantsmooth<-function(intensities, smooth.lambda=2, tau=0.5, ridge.kappa=0, smooth.na=TRUE) {
  m<-length(intensities)
  nas<-is.na(intensities)
  if (sum(nas)<m) {
    # fill missing values in with valid value
    intensities[nas]<-0
    E<-diag(m)
    Dif<-diff(E)
    # use zero weights for original missing values
    B<-rbind(diag(as.numeric(!nas)),smooth.lambda * Dif)
    ystar = c(intensities, rep(0, m - 1))
    if (ridge.kappa > 0) {
      B <- rbind(B, ridge.kappa * E)
      ystar <- c(ystar, rep(0, m))
    }
    myrq = try(rq.fit(B, ystar, tau=tau, method = "fn"),TRUE)
    if (class(myrq)!="try-error") {
      res<-myrq$coeff
      if (!smooth.na) res[nas]<-NA
      res
    }
    else {
      if (ridge.kappa==0) {
        warning("Problem with fit, repeated with ridge.kappa=(0.001*smooth.lambda)")
        .quantsmooth(intensities,smooth.lambda=smooth.lambda,tau=tau,ridge.kappa=smooth.lambda*0.001,smooth.na=smooth.na)
      }
      else {
        myrq  #Show error
      }
    }
  } else {
    warning("data is all NA, result is NA")
    rep(NA,m)
  }
}

quantsmooth.seg  <-  function(y, x = 1:length(y), lambda = 2, tau = 0.5,
                              kappa = 0, nb = length(x)) {
  # Quantile smoothing with smaller basis
  # Basis has nb segments
  # 
  # Paul Eilers, 2007
  # Based on function .quantsmooth() in package 'quantsmooth' (Oosting et al.)
                            
  # Remove NAs before computation
  nas = is.na(y)
  m = length(y)
  if (m == sum(nas)) {
      warning("data is all NA, result is NA")
      return(NA * y)
  }

  # Construct regression basis
  dx = (1+1e-7) * (max(x) - min(x)) / nb
  ix =  floor((x - min(x)) / dx) + 1
  if (nb == m) ix = 1:m
  B = outer(ix, 1:nb, '==')
  
  # Zero weights for NAs
  B[nas,] = 0
  y[nas] = 0

  # Augment data with penalty stuff
  E = diag(nb)
  D = diff(E)
  B = rbind(B, lambda * D)
  ystar = c(y, rep(0, nb - 1))
  if (kappa > 0) {
    B  =  rbind(B, kappa * E)
    ystar  =  c(ystar, rep(0, nb))
  }
  
  # Try quantile regression
  myrq = try(rq.fit(B, ystar, tau = tau, method = "fn"), TRUE)
  
  # If result is OK, return smooth result
  if (class(myrq) != "try-error") {
      a = myrq$coeff
      z = a[ix]
      return(z)
  }
  
  # If failure, even with kappa >0, return warning
  if (kappa > 0) return(myrq)
  
  # If failure and kappa >0, retry wih reasonable kappa
  warning("Problem with fit, repeated with kappa = 0.001 * lambda")
  z = quantsmooth.seg(y, x, lambda = lambda, tau = tau,
                       kappa = lambda * 0.001, nb = nb)
  return(z)
}


quantsmooth<-function(intensities, smooth.lambda=2, tau=0.5, ridge.kappa=0, smooth.na=TRUE, segment) {
  # if segment is set then the sequence is smoothed with overlapping segments
  # The algorhithm has steeply increasing memory needs for longer sequences
  m<-length(intensities)
  if (missing(segment)) segment<-m
  step.size<-segment %/% 2
  response<-vector(mode="numeric",length=m)
  response[1:min(m,segment)]<-.quantsmooth(intensities[1:min(m,segment)],smooth.lambda,tau,ridge.kappa,smooth.na)
  i.s<-1+step.size
  ol<-segment-step.size
  while ((i.s+step.size) < m) {
    i.e<-min(m,i.s+segment-1)
    tmp.resp<-.quantsmooth(intensities[i.s:i.e],smooth.lambda,tau,ridge.kappa,smooth.na)
    if (ol>0) {
      # set diagonal tapering on overlapping sequencing to prevent abrupt changes on start and end of overlap
      portion<-1:ol / (ol+1)
      response[i.s:(i.s+ol-1)] <- (response[i.s:(i.s+ol-1)]*(1-portion)) + (tmp.resp[1:ol] * portion)
      
      response[(i.s+ol):i.e]<-tmp.resp[(ol+1):length(tmp.resp)]
    } else {
      response[i.s:i.e]<-tmp.resp
    }
    i.s<-i.s+step.size 
  }
  response
}

quantsmooth.cv<-function(intensities, smooth.lambda=2, ridge.kappa=0) {
  m<-length(intensities)
  nas<-is.na(intensities)
  if (sum(nas)<m) {
    # fill missing values in with valid value
    intensities[nas]<-0

    E<-diag(m)
    Dif<-diff(E)
    ystar = as.vector(c(intensities, rep(0, m - 1)))
    weight.odd<-rep(c(1,0),length.out=m)*as.numeric(!nas)
    weight.even<-rep(c(0,1),length.out=m)*as.numeric(!nas)
  
    E.odd<-diag(weight.odd)
    B.odd<-rbind(E.odd,smooth.lambda * Dif)
    E.even<-diag(weight.even)
    B.even<-rbind(E.even,smooth.lambda * Dif)
    if (ridge.kappa > 0) {
      B.odd <- rbind(B.odd, ridge.kappa * E)
      B.even<- rbind(B.even, ridge.kappa * E)
      ystar <- c(ystar, rep(0, m))
    }
  
    myrq.odd = try(rq.fit(B.odd, ystar, method = "fn"),FALSE)
    if (class(myrq.odd)=="try-error") {
      warning("error in fit, result is NA")
      NA
    }
    else {  
      myrq.even = try(rq.fit(B.even, ystar, method = "fn"),FALSE)
      if (class(myrq.even)=="try-error") {
        warning("error in fit, result is NA")
        NA
      }
      else {  
        resid.odd<-intensities-myrq.odd$coefficients
        resid.even<-intensities-myrq.even$coefficients
        #sum of squares van interpolated values
        sum(resid.odd * resid.odd * weight.even,na.rm = TRUE) + sum(resid.even * resid.even * weight.odd,na.rm = TRUE)
      }
    }
  } else {
    warning("data is all NA, result is NA")
    NA
  }  
}
#
getLambdaMin<-function(intensities, lambdas, ...) {
  lambda.res<-rep(NA,length(lambdas))
  for (lambda in 1:length(lambdas)) lambda.res[lambda]<-quantsmooth.cv(intensities,lambdas[lambda],...)
  lambdas[which.min(lambda.res)]
}
#
plotSmoothed<-function(intensities, position, ylim=NULL, ylab="intensity", xlab="position", normalized.to=NULL, grid=NULL, smooth.lambda=2, interval=0.5, plotnew=TRUE, cols, cex.pts=0.6, ...) {
  # plot smoothed data
  # median line is drawn continuous
  # quantile intervals are plotted symmetrical around median ie interval 0.5 plots 0.25 and 0.75 quantiles
  # if intensities contains more than 1 column, the columns are drawn separately
  # position is single vector
  if(is.null(ylim)) ylim<-c(min(intensities,na.rm=TRUE),max(intensities,na.rm=TRUE))
  if(plotnew)plot(c(min(position),max(position)),ylim,ylab=ylab,xlab=xlab,type="n",...)
  if (!is.null(grid)) abline(v=grid,lty=2)
  if (!is.null(normalized.to)) abline(h=normalized.to)
  intensities<-as.matrix(intensities) # make sure it works if only a vector is supplied
	if(missing(cols)) cols<-1:ncol(intensities)+1
  
	idx<-order(position)
  position<-position[idx]
  intensities<-intensities[idx,,drop=FALSE]
  
  for (sample in 1:ncol(intensities)) {
	  if (cex.pts>0) points(position,intensities[,sample],col=cols[sample],pch=20,cex=cex.pts)
    if (sum(!is.na(intensities[,sample]))>10) {
      lines(position, quantsmooth(intensities[,sample],smooth.lambda,segment=150), col=cols[sample], lwd=2)
      if (length(interval)>0) {
        for (i in 1:length(interval)) {
          lines(position, quantsmooth(intensities[,sample],smooth.lambda,tau=0.5-(interval[i]/2),segment=150), col=cols[sample], lty=1+i)
          lines(position, quantsmooth(intensities[,sample],smooth.lambda,tau=0.5+(interval[i]/2),segment=150), col=cols[sample], lty=1+i)
        }
      }
    }
  }
}

getChangedIdx<-function(changed, up) {
  if (sum(changed,na.rm=TRUE)>0) {
    changed[is.na(changed)]<-FALSE
    crossing<-xor(c(FALSE,changed),c(changed,FALSE))
    position<-seq(1,length.out=length(crossing))[crossing]
    startidx<-seq(1,by=2,length.out=length(position) / 2) # odd indexes
    startpos<-position[startidx]
    endpos<-position[startidx+1]-1
    data.frame(up=up,start=startpos,end=endpos)
  } else NULL
}

getChangedRegions<-function(intensities, positions, normalized.to=1, interval, threshold, minlength=2, ...) {
  # determine regions with changes after smoothing
  # normalized.to: value to compare with
  # smooth.lambda: smoothing parameter
  # interval     : changes are defined by these smoothed boundaries crossing normalized.to
  # treshold     : changes are defined by croosing of signal outside of normalized.to + or - treshold 
  #                (only one of treshold or interval can be defined)
  # minlength    :  minimum length of a change to be listed
  #
  # value        : dataframe 3 columns up, start, end
	if (missing(positions)) positions<-1:length(intensities)
	if (!is.null(match.call()$tau)) stop("tau is set by the function")
	if (length(positions)!=length(intensities)) stop("Length of positions argument should be equal to length of intensities argument")
  if (!missing(interval)) {
    res<-rbind(getChangedIdx(quantsmooth(intensities,tau=0.5-(interval/2),...) > normalized.to,TRUE),
          getChangedIdx(quantsmooth(intensities,tau=0.5+(interval/2),...) < normalized.to,FALSE))
  } else if (!missing(threshold)) {
    smoothed<-quantsmooth(intensities,tau=0.5,...)
    res<-rbind(getChangedIdx(smoothed > (normalized.to+threshold),TRUE),
          getChangedIdx(smoothed < (normalized.to-threshold),FALSE))
  } else stop("Either treshold or interval should be defined")
  if (!is.null(res)) {
	  res[,"start"]<-positions[res[,"start"]]
	  res[,"end"]<-positions[res[,"end"]]
	} 
	res 
}

