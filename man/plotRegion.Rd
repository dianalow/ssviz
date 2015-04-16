\name{plotRegion}
\alias{plotRegion}
\title{plotRegion
}
\description{
Plots the read density given a chromosome region.
}
\usage{
plotRegion(bamlist, region, howsmooth = 2, ncounts = NULL, samplenames = NULL)
}
\arguments{
  \item{bamlist}{An object of type list, giving a list of bam files. If you only have 1 file, use \code{list(bam_file)}
}
  \item{region}{An object of type character defining the region to plot. Eg. \code{chr1:1000-2000}
}
  \item{howsmooth}{Numeric value controlling smoothness of the plot.
}

  \item{ncounts}{Total number of reads for plot normalization.
}
  \item{samplenames}{Sample names
}
}
\value{Returns the x and y components of the region's reads and plots the density.
}
\author{
Diana H.P. Low
}
\examples{
data(ssviz)
region<-'chr1:3015526-3080526'
plotRegion(list(ctrlbam), region=region)
}
