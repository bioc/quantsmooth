\name{numericCHR}
\alias{numericCHR}
\alias{characterCHR}
\title{ Conversion of chromosome IDs between numeric and character }
\description{
  The function converts chromosomal ids to their numeric form, and the sex 
  chromosomes to values between 98 and 100.
  This simplifies sorting on chromosome ID
}
\usage{
numericCHR(CHR)
characterCHR(CHR)
}
\arguments{
  \item{CHR}{character/numeric vector for both functions the mode of the input 
              is not forced. For numericCHR strings "X","Y" and "XY" are converted 
              to 98,99 and 100 respectively.  }
}
\value{
  \code{numericCHR} returns a numeric vector of same length as \code{CHR}
  \code{characterCHR} returns a character vector of same length as \code{CHR}
}
\author{ Jan Oosting }
\examples{
}
\keyword{manip}
