\name{pingpong}
\alias{pingpong}
\title{
pingpong
}
\description{
piRNA ping-pong analysis of complementary sequences
}
\usage{
pingpong(bam_file)
}
\arguments{
  \item{bam_file}{An object of class data.frame or DataFrame
}
}
\details{The ping-pong mechanism is a proposed method for the amplification of primary piRNAs, which leads to the production of new primary piRNAs from their precursor transcripts, which eventually amplifies the pool of both primary and secondary piRNAs. This positive feedback loop is a secondary biogenesis mechanism that requires complementary transcripts to a pre-existing pool of piRNAs.
}
\value{This function returns a data.frame object with frequency of overlapping complementary piRNAs.
}
\references{
Brennecke J. et al. Cell 128, 1089-1103, March 23, 2007
}
\author{
Diana H.P. Low
}

\examples{
data(ssviz)
pp<-pingpong(pctrlbam)
}