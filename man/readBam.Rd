\name{readBam}
\alias{readBam}
\title{
readBam
}
\description{
Reads a bam file through RSamtools, and converts it into a data frame of class \code{DataFrame}
}
\usage{
readBam(file_name, tags = character(0))
}
\arguments{
  \item{file_name}{Character string of bam file location
}
  \item{tags}{Bam tags to import into the data frame. By default it only takes the standard values if none are given.
}
}
\details{
This function formalizes what had been described in the RSamtools documentation and makes it easier to compute the downstream functions in this package.
}
\value{Returns the bam file contents in a readable dataframe format.
}
\references{
RSamtools package
}
\author{
Diana H.P. Low
}

\examples{
bam.files <- dir(system.file("extdata", package = "ssviz"), full = TRUE, patt = "bam$")
ctrlbam <- readBam(bam.files[1])
}
