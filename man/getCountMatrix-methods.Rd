\name{getCountMatrix-methods}
\docType{methods}
\alias{getCountMatrix-methods}
\alias{getCountMatrix,DataFrame-method}
\title{getCountMatrix}
\description{
returns the bam data.frame with an additional column counts. Only relevant if the fasta file used for mapping input was previously collapsed via fastx_toolkit to return a fasta read name in the format of \code{readnumber-totalcounts}
}
\section{Methods}{
\describe{

\item{\code{signature(object="DataFrame")}}{Returns and object of class data.frame having the values from the original bam file with an additional 'count' column.
}
}}