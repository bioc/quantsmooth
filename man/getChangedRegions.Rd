\name{getChangedRegions}
\alias{getChangedRegions}

\title{getChangedRegions}

\description{
   retrieve regions of interest in a vector of intensities using quantile smoothing
}

\usage{
  getChangedRegions(intensities, positions, normalized.to=1, smooth.lambda=2, interval, threshold, minlength=2)
}

\arguments{
   \item{intensities}{numeric vector}
   \item{positions}{numeric vector of the same length as intensities. If this argument is not given the results contain the indexes of the \code{intensities} vector, else the values in \code{positions} are used. Both vectors are sorted in the order of \code{positions}. }
   \item{normalized.to}{numeric, reference value. Changes are compared to this value}
   \item{smooth.lambda}{smoothing parameter}
   \item{interval}{numeric [0,1], bandwidth around reference. If the smoothed line at the higher quantile drops below the \code{normalized.to} value, a deleted region is recognized, and vice versa.}
   \item{threshold}{numeric, if the median smoothed value drops below \code{normalized.to - threshold}, or above \code{normalized.to + threshold} a changed region is called}
   \item{minlength}{integer, not used currently}
}

\details{
  This function uses \code{quantsmooth} to detect regions in the genome that are abnormal.
  If \code{interval} is set then a smoothed line is calculated for \code{tau = 0.5 - interval/2}, and a region is determined as upregulated if this line is above the reference. Down regulation is determined when the smoothed line for \code{tau = 0.5 + interval/2} is below the reference value.
  If \code{threshold} is set then a smoothed line is calculated for \code{tau = 0.5} and up- or down regulation are determined when this line is outside the range \code{[normalized.t - threshold:normalized.to + threshold]}
}
\value{
  A data.frame with 3 colums is returned. Each row contains a region with columns \code{up}, \code{start} and \code{end}. \code{start} and \code{end} indicate positions in the vector of the first and last position that were up- or downregulated
}
\author{Jan Oosting}

\examples{
  data(chr14)
  getChangedRegions(ill.cn[,1],ill.pos,normalized.to=2,interval=0.5)
}
\keyword{attribute}

