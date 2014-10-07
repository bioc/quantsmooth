\name{grid.chromosome}
\alias{grid.chromosome}
\title{ Draw a chromosome using the grid package }
\description{
  A chromosme is drawn including the cytobands
}
\usage{
grid.chromosome(chrom, side = 1, units = "hg19", chrom.width = 0.5, length.out, 
             bands = "major", legend = c("chrom", "band", "none"), cex.leg = 0.7, bleach = 0, ...)
}
\arguments{
  \item{chrom}{numeric or character, id of chromosome to plot}
  \item{side}{numeric [1:4], side of rectangle to draw, 4 sides, side 2 and 4 are vertical}
  \item{units}{character or data.frame, type of units for genomic data, or a dataframe with UCSC cytoband data, see \code{\link{lengthChromosome}}}
  \item{chrom.width}{numeric [0,1], The width relative to the width (sides 2 and 4) or height(sides 1 and 3) of the viewport}
  \item{length.out}{numeric, size of native units of viewport}
  \item{bands}{character, draw either major or minor bands}
  \item{legend}{character, type of legend}
  \item{cex.leg}{numeric, relative size of legend text}
  \item{bleach}{numeric [0,1], proportion by which to bleach the chromosome}
  \item{\dots}{arguments for viewport(), especially x,y, width, and height}
}
\details{
   The chromosome is drawn within a rectangle defined by x, y, width, and height, which is pushed as a viewport.
The legend is drawn within the same rectangle in the space left over by chrom.width.    
}
\references{ lodplot package }
\value{
This function is executed for its side effects
}

\author{David L Duffy ,Jan Oosting}

\seealso{\code{\link{paintCytobands}}}
\examples{
    
    grid.newpage()
    grid.chromosome(1,units="bases",height=0.15)

}
\keyword{aplot}