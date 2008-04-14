\name{helperFunctions}
\alias{getChangedIdx}

\title{arraysmooth Helper functions}

\description{
   quantsmooth Helper functions
}

\usage{
  getChangedIdx(changed,up)
}

\arguments{
   \item{changed}{}
   \item{up}{}
}

\details{
  Should not be called directly
}
\value{
  A data.frame with 3 colums is returned. Each row contains a region with columns \code{up}, \code{start} and \code{end}. \code{start} and \code{end} indicate positions in the vector of the first and last position that were up- or downregulated
}
\author{Jan Oosting}

\keyword{internal}
