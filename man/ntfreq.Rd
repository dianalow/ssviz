\name{ntfreq}
\alias{ntfreq}
\title{
ntfreq
}
\description{
Calculates nucleotide frequency of reads in bam file
}
\usage{
ntfreq(bam_file, ntlength, toRNA = TRUE, count_type = "total")
}
\arguments{
  \item{bam_file}{An object of class data.frame or DataFrame
}
  \item{ntlength}{ An integer specifying the length of the sequence to quantify
}
  \item{toRNA}{A logical value on whether to translate the DNA sequence to RNA
}
  \item{count_type}{A character string on how to count the nucleotides. Can be either "total" or "unique". If \code{total} is selected, the function will look for the \code{count}column and multiply the reads by its number of occurence when calculating the frequency.
}
}
\value{
Returns a data.frame of the frequency of nucleotides (either A/C/G/T or A/C/G/U) at each position up to the specified \code{ntlength}
}
\author{
Diana H.P. Low
}

\examples{
data(ssviz)
freq<-ntfreq(pctrlbam,ntlength=10)
}
