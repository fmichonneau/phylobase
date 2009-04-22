\name{hasSingle}
\alias{hasSingle}
\alias{hasPoly}
\alias{hasRetic}
\title{Test trees for polytomies, inline nodes,
or reticulation}
\description{checks to see whether trees have
  (structural) polytomies, inline nodes (i.e.,
  nodes with a single descendant), or reticulation
  (i.e., nodes with more than one ancestor)
}
\usage{
hasSingle(object)
hasPoly(object)
hasRetic(object)
}
\arguments{
  \item{object}{an object inheriting from class \code{phylo4}}
}
\value{
  Logical value
}
\author{Ben Bolker}
\note{
  Some algorithms are unhappy with structural polytomies (i.e., >2
  descendants from a node), with single-descendant nodes,
  or with reticulation; these
  functions check those properties.
  We haven't bothered to check for zero branch lengths:
  the consensus is that it doesn't come up much,
  and that it's simple enough to
  test \code{any(edgeLength(x) == 0)} in these cases.
  (Single-descendant nodes are used e.g. in OUCH, or in
  other cases to represent events occurring along a branch.)
}
\examples{
library(ape)
tree.owls.bis <- read.tree(text = "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);")
owls4 <- as(tree.owls.bis, "phylo4")
hasPoly(owls4)
hasSingle(owls4)
}
\keyword{misc}
