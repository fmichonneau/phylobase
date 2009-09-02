\name{reorder-methods}
\docType{methods}
\alias{reorder-methods}
\alias{reorder,phylo-method}
\alias{reorder,phylo4-method}
\alias{reorder,phylo4d-method}
\title{reordering trees within phylobase objects}
\description{
 Methods for reordering trees into various traversal orders
}
\section{Methods}{
  \describe{
    \item{x = "phylo"}{reorders a \code{phylo} object}
    \item{x = "phylo4"}{reorders a \linkS4class{phylo4} object}
    \item{x = "phylo4d"}{reorders a \linkS4class{phylo4d} object}
  }
}

\usage{
\S4method{reorder}{phylo}(x, order = "cladewise")
\S4method{reorder}{phylo4}(x, order = c("preorder", "postorder"))
}

\arguments{
  \item{x}{a \code{phylo4} or \code{phylo4d} object}
  \item{order}{The desired traversal order; currently only 'preorder'
    and 'postorder' are allowed for \code{phylo4} and \code{phylo4d}
    objects, whereas only 'cladewise' and 'pruningwise' are allowed for
    \code{phylo} objects}
}
\details{
  The \code{reorder} method takes a \code{phylo4} or \code{phylo4d} tree
  and orders the edge matrix (i.e. \code{edges(x)}) in the requested
  traversal order. Currently only two orderings are permitted, and both
  require rooted trees. In "postorder", a node's descendants come before
  that node, thus the root, which is ancestral to all nodes, comes last.
  In "preorder", a node is visited before its descendants, thus the root
  comes first.

  A method is also defined that takes an \code{ape phylo} object.  This
  also takes an order argument, however, 'pruningwise' and 'cladewise'
  are the only acceptable parameters. This is because this method
  actually uses the \code{ape reorder()} command to complete the
  ordering.
}
\value{
  A \code{phylo4} or \code{phylo4d} object with the edge, label, length
  and data slots ordered as \code{order}, which is itself recorded in
  the order slot.
}
\author{Peter Cowan, Jim Regetz}
\note{
  The "preorder" parameter corresponds to "cladewise" in the \code{ape}
  package, and "postorder" corresponds (almost but close enough?) to
  "pruningwise".  
  
  See \url{http://ape.mpl.ird.fr/misc/FormatTreeR_28July2008.pdf}
}
\seealso{\code{\link[ape]{reorder.phylo}} in the \code{ape} package. 
 \code{\link{ancestors}} \code{\link{ancestor}} \code{\link{siblings}}
 \code{\link{children}} \code{\link{descendants}}
}
\examples{
phy <- phylo4(rtree(5))
edges(reorder(phy, "preorder"))
edges(reorder(phy, "postorder"))
}
\keyword{methods}
