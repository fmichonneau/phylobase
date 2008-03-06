\name{hasSingles}
\alias{hasSingles}
\alias{hasPoly}
\title{Test trees for polytomies or inline nodes}
\description{checks to see whether trees have
  (structural) polytomies or inline nodes (i.e.,
  nodes with a single descendant)
}
\usage{
hasSingles(object)
hasPoly(object)
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
  descendants from a node) or with single-descendant nodes; these
  functions check those properties.
  We haven't bothered to check for zero branch lengths:
  the consensus is that it doesn't come up much,
  and that it's simple enough to
  test \code{any(edgeLength(x)==0)} in these cases.
  (Single-descendant nodes are used e.g. in OUCH, or in
  other cases to represent events occurring along a branch.)
}
\examples{
library(ape)
example(read.tree)
owls4 = as(tree.owls.bis,"phylo4")
hasPoly(owls4)
hasSingles(owls4)
}
\keyword{misc}
