\name{quantsmooth.seg}
\alias{quantsmooth.seg}

\title{quantsmooth.seg}

\description{
   segmented Quantile smoothing of array data
}

\usage{
  quantsmooth.seg(y, x = 1:length(y), lambda = 2, tau = 0.5,kappa = 0, nb = length(x))
}

\arguments{
   \item{y}{numeric vector}
   \item{x}{numeric vector of same length as \code{y}. Position of values}
   \item{lambda}{numeric}
   \item{tau}{numeric [0..1], the quantile desired; see \code{\link[quantreg]{rq.fit}}}
   \item{kappa}{fudge parameter; see details}
   \item{nb}{integer, basis}
}

\details{

}
\value{
  This function returns a vector of the same length as \code{y}
}
\author{Jan Oosting}

\examples{
	data(chr14)
	plot(quantsmooth.seg(bac.cn[,1],lambda=2.8,nb=50),type="l")
}
\keyword{smooth}
