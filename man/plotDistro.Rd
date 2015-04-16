\name{plotDistro}
\alias{plotDistro}
\title{
plotDistro
}
\description{
Plots distribution of reads in the bam file based on length, direction (strand) or location (rname)
}
\usage{
plotDistro(bamlist, type = "qwidth", samplenames = NULL, unique = FALSE, ncounts = NULL, norm = FALSE, yname = "Counts per million")
}

\arguments{
  \item{bamlist}{ An object of type list, giving a list of bam files. If you only have 1 file, use \code{list(bam_file)}
}
  \item{type}{ An object of type character. Can be \code{qwidth,rname} or \code{strand}. In theory, any column property existing in the bam file can be used, but these 3 would be most meaningful.
}
  \item{samplenames}{Labels for the plot.
}
  \item{unique}{ Logical value to use unique reads (\code{TRUE}) or all reads (\code{FALSE})
}
  \item{ncounts}{Number of total counts in the bam file, used if \code{unique} is set to \code{FALSE}.
}
  \item{norm}{ Logical value to determine if plot will be normalised.
}
  \item{yname}{y axis label.
}
}
\author{
Diana H.P. Low
}
\examples{
data(ssviz)
plotDistro(list(ctrlbam))
}
