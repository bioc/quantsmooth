\name{Chromosome14}
\alias{chr14}
\alias{affy.cn}
\alias{affy.pos}
\alias{ill.cn}
\alias{ill.pos}
\alias{bac.cn}
\alias{bac.pos}

\docType{data}

\title{Example data from several quantitative genomic methods}

\description{ A collection of arrays that contains data of chromosome 14 of 3 colorectal tumors.
The first tumor shows 1 region of loss, the second tumor shows no abberation, while the third 
tumor shows loss of 1 copy of the chromosome.
  \describe{
   \item{affy.cn}{ Copy number values of 358 probes from Affymetrix 10K genechip. Data was obtained from
    DChip}
   \item{affy.pos}{ corresponding probe positions}
   \item{bac.cn}{ Copy number values of 112 probes from a 1 mb spaced BAC array-CGH }
   \item{bac.pos}{ corresponding probe positions}
   \item{ill.cn}{ Copy number values of 207 probes from Illumina GoldenGate Linkage IV data}
   \item{ill.pos}{ corresponding probe positions}
  }
}

\usage{
   data(chr14)
}

\format{Matrices of copy number values and vectors of chromosomal probe positions}

\author{Jan Oosting}

\keyword{datasets}
