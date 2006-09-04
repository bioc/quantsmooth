\name{quantsmooth.cv}
\alias{quantsmooth.cv}

\title{quantsmooth.cv}

\description{
   Cross validation of smoothing parameters
}

\usage{
  quantsmooth.cv(intensities,smooth.lambda=2, ridge.kappa=0)
}

\arguments{
   \item{intensities}{numeric vector}
   \item{smooth.lambda}{numeric; see \code{\link{quantsmooth}}}
   \item{ridge.kappa}{fudge parameter; see \code{\link{quantsmooth}}}
}

\details{
  Cross validation is performed by calculating the fit from the even indices on the odd indices and vice versa.
}
\value{
  This function returns the sum of squared differences or \code{NA} if the fitting function gave an error
}
\author{Jan Oosting}

\seealso{
  \code{\link{getLambdaMin}}
}

\examples{
	data(chr14)
	# A low value is indicative of a better fit to the data
	quantsmooth.cv(bac.cn[,1],1)
	quantsmooth.cv(bac.cn[,1],2.8)
}
\keyword{htest}
\keyword{smooth}
