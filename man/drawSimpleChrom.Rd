\name{drawSimpleChrom}
\alias{drawSimpleChrom}
\title{ Draw chromosome-like icons }
\description{
  This function paints chromosomal icons on an existing plot
}
\usage{
drawSimpleChrom(x, y, len = 3, width = 1, fill, col, orientation = c("h", "v"), centromere.size = 3 * width)
}
\arguments{
  \item{x}{start x-position}
  \item{y}{start y-position}
  \item{len}{total length of the chromosome}
  \item{width}{width of the chromosome}
  \item{fill}{character, \{"a","p","q","q[1-3]","p[1-3]" \}. Events to a chromosome can be depicted by coloring
    "a"ll of the chromosome, the complete p or q-arm, or a subsegment of the arms}
  \item{col}{color(s) of \code{fill}}
  \item{orientation}{either "h"orizontal or "v"ertical}
  \item{centromere.size}{The centromere is drawn with a point, and may need a size adjustment}
}
\value{
  This function is executed for its side effects
}
\author{ Jan Oosting }
\examples{
}
\keyword{aplot}
