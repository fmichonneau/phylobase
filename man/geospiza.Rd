\name{geospiza}
\alias{geospiza}
\alias{geospiza_raw}
\docType{data}
\title{Data from Darwin's finches}
\description{
  Phylogenetic tree and morphological data for Darwin's finches,
in different formats
}
\usage{
   data(geospiza)
   data(geospiza_raw)
}
\format{
  \code{geospiza} is a \code{phylo4d} object; 
\code{geospiza_raw} is a list containing \code{tree}, a \code{phylo}
object (the tree), \code{data}, and a data frame with the data
(for showing examples of how to merge tree and data)
}
\note{
  Stolen from Luke Harmon's Geiger package,
  to avoid unnecessary dependencies
}
\source{
  Dolph Schluter via Luke Harmon
}
\examples{
data(geospiza)
plot(geospiza,adj=.4,cex.lab=.7,cex.sym=.8)
}
\keyword{datasets}
