\name{reorder-methods}
\docType{methods}
\alias{reorder-methods}
\alias{reorder,phylo-method}
\alias{reorder,phylo4-method}
\alias{reorder,phylo4d-method}
\alias{orderIndex}
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
\S4method{reorder}{phylo}(x, ...)
\S4method{reorder}{phylo4}(x, ...)
\S4method{reorder}{phylo4d}(x, ...)
orderIndex(phy, order = c('preorder','postorder'))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a \code{phylo4} or \code{phylo4d} object}
  \item{x}{a \code{phylo4} or \code{phylo4d} object}
  \item{order}{The desired traversal order currently only 'preorder' and 'postorder' are allowed}
  \item{\dots}{ A required order agrument (see above) is the only additional argument accepted. }
}
\details{
The \code{reorder} method takes \code{phylo4} or \code{phylo4d} tree and orders the edge matrix (e.g. \code{edges(phy)}) in the requested traversal order.  Currently only two orderings are permitted and both require rooted trees.  In "postorder" a nodes descendants come before that node, thus the root, being ancestral to all nodes, comes last.  "Preorder" is the converse, a node is visited before it's descendants, thus the root comes first.

A method is also defined that takes an \code{ape phylo} object.  This also takes an order argument, however, 'pruningwise' and 'cladewise' are the only acceptable parameters.  This is because this method actually uses the \code{ape reorder()} command to complete the ordering.  
}
\value{
A \code{phylo4} or \code{phylo4d} object with the edge, label, lenght and data slots ordered as \code{order}.
}
\author{Peter Cowan}
\note{
  The "preorder" parameter corresponds to "cladewise" in the \code{ape} package, and "postorder" corresponds (almost but close enough?) to "pruningwise".  
  
  See \url{http://ape.mpl.ird.fr/misc/FormatTreeR_28July2008.pdf}
}
\seealso{\code{\link[ape]{reorder.phylo}} in the \code{ape} package. 
 \code{\link{ancestors}} \code{\link{ancestor}} \code{\link{siblings}}
 \code{\link{children}} \code{\link{descendants}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
