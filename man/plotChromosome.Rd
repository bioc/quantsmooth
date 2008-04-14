\name{plotChromosome}
\alias{plotChromosome}
\title{Wrapper for plotSmoothed}
\description{
  This function is a wrapper for plotSmoothed, to make data subsetting easier
}
\usage{
plotChromosome(gendata, chrompos, chromosome, dataselection = NULL, ylim = NULL, normalized.to = NULL, grid = NULL, smooth.lambda = 2, interval = 0.5, ...)
}
\arguments{
  \item{gendata}{numeric matrix or data.frame}
  \item{chrompos}{chrompos object with same numer of rows as gendata}
  \item{chromosome}{numeric, chromosme to show}
  \item{dataselection}{optional, subset of samples/columns in gendata}
  \item{ylim}{limits for plot}
  \item{normalized.to}{y-value(s) for line}
  \item{grid}{x-value(s) for line}
  \item{smooth.lambda}{smoothing parameter, see \code{\link{quantsmooth}}}
  \item{interval}{position of extra lines besides median, see \code{\link{plotSmoothed}}}
  \item{\dots}{extra arguments for \code{\link{plotSmoothed}}}
}
\value{
  The function is used for its side effects
}
\author{ Jan Oosting }
\seealso{ \code{\link{plotSmoothed}}, \code{\link{quantsmooth}} }
\keyword{hplot}
