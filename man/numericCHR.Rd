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
numericCHR(CHR, prefix="chr")
characterCHR(CHR, prefix="")
}
\arguments{
  \item{CHR}{character/numeric vector for both functions the mode of the input 
              is not forced. For numericCHR strings "X","Y" and "XY" are converted 
              to 98,99 and 100 respectively.  }
  \item{prefix}{character, string is excluded from (\code{numericCHR}) or prepended to (\code{characterCHR}) all items of the output}
}
\value{
  \code{numericCHR} returns a numeric vector of same length as \code{CHR}
  \code{characterCHR} returns a character vector of same length as \code{CHR}
}
\author{ Jan Oosting }
\examples{
   chroms<-c("3","2","8","X","7","Y","5","1","9","10","11","12","4","6")
   sort(chroms)
   sort(numericCHR(chroms))
   characterCHR(sort(numericCHR(chroms)),prefix="chr")
}
\keyword{manip}
