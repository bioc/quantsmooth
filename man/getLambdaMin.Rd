\name{getLambdaMin}
\alias{getLambdaMin}

\title{getLambdaMin}

\description{
   Test a set of smoothing parameters to find best fit to data
}

\usage{
  getLambdaMin(intensities,lambdas,\dots)
}

\arguments{
   \item{intensities}{numeric vector}
   \item{lambdas}{numeric vector; see \code{\link{quantsmooth}}}
   \item{\dots}{extra parameters for \code{quantsmooth.cv}; currently only \code{ridge.kappa}}
}

\details{
  Cross validation is performed using a set of lambda values in order to find the lambda value that shows the best fit to the data.
}
\value{
  This function returns the lambda value that has the lowest cross validation value on this dataset
}
\author{Jan Oosting}

\seealso{
  \code{\link{quantsmooth.cv}}
}

\examples{

}
\keyword{smooth}
