\name{getCountMatrix}
\alias{getCountMatrix}
\title{
getCountMatrix
}
\description{
returns the bam data.frame with an additional column counts. Only relevant if the fasta file used for mapping input was previously collapsed via fastx_toolkit to return a fasta read name in the format of \code{readnumber-totalcounts}
}
\usage{
getCountMatrix(bam_file)
}
\arguments{
  \item{bam_file}{An object of class DataFrame (from IRanges). Can be generated from readBam.
}
}
\value{An object of class data.frame having the values from the original bam file with an additional 'count' column.
}
\author{ Diana H.P. Low
}

\seealso{
\code{\link{readBam}}
}
\examples{
data(ssviz)
getCountMatrix(ctrlbam)
}
