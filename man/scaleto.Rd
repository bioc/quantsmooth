\name{scaleto}
\alias{scaleto}
\title{ Scales data within a range to a new range}
\description{
  This function scales data to a new range while enforcing the boundaries. This 
  can be helpful in preventing overlap between chromosomal plots that display
  multiple chromosomes in the same plot
}
\usage{
scaleto(x, fromlimits = c(0, 50), tolimits = c(0.5, -0.5), adjust = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric}
  \item{fromlimits}{numeric vector with length 2, original range of data}
  \item{tolimits}{numeric vector with length 2, target range of data}
  \item{adjust}{logical, if \code{TRUE} then the target values are clipped to
                the target range}
}
\value{
  numeric of same size as \code{x}
}
\author{ Jan Oosting }
\examples{
}
\keyword{manip}
