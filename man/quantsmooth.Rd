\name{quantsmooth}
\alias{quantsmooth}

\title{quantsmooth}

\description{
   Quantile smoothing of array data
}

\usage{
  quantsmooth(intensities,smooth.lambda=2, tau=0.5, ridge.kappa=0,smooth.na=TRUE,segment)
}

\arguments{
   \item{intensities}{numeric vector}
   \item{smooth.lambda}{numeric}
   \item{tau}{numeric [0..1], the quantile desired; see \code{\link[quantreg]{rq.fit}}}
   \item{ridge.kappa}{fudge parameter; see details}
   \item{smooth.na}{logical; handling of NA}
   \item{segment}{integer, length of overlapping segments}
}

\details{

}
\value{
  This function returns a vector of the same length as \code{intensities}, or a matrix
  if the length of \code{tau} is greater than 1.
}
\author{Jan Oosting}

\examples{
	data(chr14)
	plot(quantsmooth(bac.cn[,1],smooth.lambda=2.8),type="l")
}
\keyword{smooth}
