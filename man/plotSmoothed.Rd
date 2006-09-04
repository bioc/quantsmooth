\name{plotSmoothed}
\alias{plotSmoothed}

\title{plotSmoothed}

\description{
  Plot a smoothed line together with the original data values
}

\usage{
  plotSmoothed(intensities, position, ylim=NULL, ylab="intensity", xlab="position", normalized.to=NULL, grid=NULL, smooth.lambda=2, interval=0.5, plotnew=TRUE, cols, \dots)
}

\arguments{
   \item{intensities}{numeric vector or matrix, data are plotted by column}
   \item{position}{numeric vector; the length should be the number of rows in intensities}
   \item{ylim}{numeric vector of length 2, limits for plot. If \code{NULL} then the minimal and maximal value in \code{intensities} is used}
   \item{ylab}{character, label for y-position}
   \item{xlab}{character, label for x-position}
   \item{normalized.to}{numeric, a line(s) is drawn at this horizontal position}
   \item{grid}{numeric, a line(s) is drawn at this vertical position}
   \item{smooth.lambda}{numeric, smoothing parameter see \code{quantsmooth}}
   \item{interval}{numeric (0..1), plotting of extra smoothed lines around median. With \code{interval} = 0.5 the 0.25 and 0.75 quartiles are plotted, with \code{interval} = 0.9 the 0.05 and 0.95 quantiles are plotted, }
   \item{plotnew}{logical, if TRUE a new plot is created, else the data are plotted into an existing plot}
   \item{cols}{color vector, colors for columns in \code{intensities} }
   \item{\dots}{extra parameters for \code{plot}}
}

\details{
  This function plots the raw data values as dots and the median smoothed values as a continuous line. If interval is supplied these are plotted as lines in different line types. More than 1 interval can be given.
  
}
\value{
  This function is used for its side effects
}
\author{Jan Oosting}

\seealso{
  \code{\link{quantsmooth}}
}

\examples{
	 data(chr14)
	 plotSmoothed(bac.cn,bac.pos,ylim=c(1,2.5),normalized.to=2,smooth.lambda=2.5)
}
\keyword{hplot}
