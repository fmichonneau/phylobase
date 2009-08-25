\name{phylo4-methods}
\docType{methods}
\alias{phylo4}
\alias{phylo4-methods}
\alias{phylo4,matrix-method}
\alias{phylo4,phylo-method}
\title{Create a phylogenetic tree}
\description{
  \code{phylo4} is a generic constructor that creates a phylogenetic tree object for use in phylobase methods. Phylobase contains functions for input of phylogenetic trees and data, manipulation of these objects including pruning and subsetting, and plotting. The phylobase package also contains translation functions to forms used in other comparative phylogenetic method packages.
}
\section{Methods}{
  \describe{
     \item{x = "matrix"}{creates a phylobase tree from a matrix of
       edges}

     \item{x = "phylo"}{creates a phylobase tree from an object of class
       \code{phylo}}
   }}
\usage{
     \S4method{phylo4}{matrix}(x, edge.length = NULL, tip.label = NULL, node.label = NULL,
        edge.label = NULL,  order="unknown", annote=list(), ...)
     \S4method{phylo4}{phylo}(x, check.node.labels = c("keep", "drop"), annote=list())
}
\arguments{
  \item{x}{a matrix of edges or an object of class \code{phylo} (see
    above)}
  \item{edge}{A numeric, two-column matrix with as many rows as branches in the phylogeny.}
  \item{edge.length}{Edge (branch) length. (Optional)}
  \item{tip.label}{A character vector of species names (names of "tip" nodes). (Optional)}
  \item{node.label}{A character vector of internal node names. (Optional)}
  \item{edge.label}{A character vector of edge (branch) names. (Optional)}
  \item{order}{character: tree ordering (allowable values are
      listed in \code{phylo4_orderings}, currently "unknown", "preorder"
      (="cladewise" in \code{ape}), and "postorder", with "cladewise"
      and "pruningwise" also allowed for compatibility with \code{ape})}
  \item{check.node.labels}{if \code{x} is of class \code{phylo}, either
    "keep" (the default) or "drop" node labels. This argument is useful
    if the \code{phylo} object has non-unique node labels.}
  \item{annote}{any additional annotation data to be passed to the new object}
  \item{\dots}{Additional arguments} %TODO a more elaborate description would be useful
}
\details{
  The minimum information necessary to create a phylobase tree object is a valid edge matrix. The edge matrix describes the topology of the phylogeny. Each row describes a branch of the phylogeny, with the (descendant) node number in column 2 and its ancestor's node number in column 1. These numbers are used internally and must be unique for each node.
  The labels designate either nodes or edges. The vector \code{node.label} names internal nodes, and together with \code{tip.label}, name all nodes in the tree. The vector \code{edge.label} names all branches in the tree. All label vectors are optional, and if they are not given, internally-generated labels will be assigned. The labels, whether user-specified or internally generated, must be unique as they are used to join species data with phylogenetic trees.
}
\author{phylobase team}
\note{Translation functions are available from many valid tree formats. See \link{coerce-methods}.  
}
\seealso{\code{\link{coerce-methods}} for translation functions. The \linkS4class{phylo4} class, the \code{\link{formatData}}
  function to check the validity of \code{phylo4} objects. See also the \code{\link{phylo4d}} constructor, and \linkS4class{phylo4d} class.}
\examples{
# a three species tree:
mytree <- phylo4(x=matrix(data=c(4,1, 4,5, 5,2, 5,3,NA,4), ncol=2, byrow=TRUE), tip.label=c("speciesA", "speciesB", "speciesC"))
mytree
plot(mytree)

# another way to specify the same tree:
mytree <- phylo4(x=cbind(c(4,4,5,5,NA), c(1,5,2,3,4)), tip.label=c("speciesA", "speciesB", "speciesC"))

# another way:
mytree <- phylo4(x=rbind(c(4,1), c(4,5), c(5,2), c(5,3), c(NA,4)), tip.label=c("speciesA", "speciesB", "speciesC"))

# with branch lengths:
mytree <- phylo4(x=rbind(c(4,1), c(4,5), c(5,2), c(5,3), c(NA,4)),
tip.label=c("speciesA", "speciesB", "speciesC"), edge.length=c(1, .2,
.8, .8, NA))
plot(mytree)
}
\keyword{classes}
