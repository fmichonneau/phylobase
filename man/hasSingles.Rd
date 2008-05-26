\name{hasSingles}
\alias{hasSingles}
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
hasSingles(object)
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
example(read.tree)
owls4 = as(tree.owls.bis, "phylo4")
hasPoly(owls4)
hasSingles(owls4)
}
\keyword{misc}
