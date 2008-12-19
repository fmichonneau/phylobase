\name{phylo4}
\alias{phylo4}
\title{ create a phylogenetic tree object }
\description{
  The \code{phylo4} function creates a phylogenetic tree object for use in phylobase methods. Phylobase contains functions for input of phylogenetic trees and data, manipulation of these objects including pruning and subsetting, and plotting. The phylobase package also contains translation functions to forms used in other comparative phylogenetic method packages.
}
\usage{
phylo4(edge, edge.length = NULL, tip.label = NULL, node.label = NULL,
edge.label = NULL, root.edge = NULL, order="unknown", ...)
}
\arguments{
  \item{edge}{A numeric, two-column matrix with as many rows as branches in the phylogeny.}
  \item{edge.length}{Edge (branch) length. (Optional)}
  \item{tip.label}{A character vector of species names (names of "tip" nodes). (Optional)}
  \item{node.label}{A character vector of internal node names. (Optional)}
  \item{edge.label}{A character vector of edge (branch) names. (Optional)}
  \item{root.edge}{Edge (branch) length from the most basal node to it's
    ancestor. (Optional)}
  \item{order}{character: tree ordering (allowable values are
      listed in \code{phylo4_ordering}, currently "unknown", "reorder"
      (="pruningwise" in \code{ape}), "postorder" (="cladewise" in \code{ape}))}
}
\details{
  The minimum information necessary to create a phylobase tree object is a valid edge matrix. The edge matrix describes the topology of the phylogeny. Each row describes a branch of the phylogeny, with the (descendant) node number in column 2 and it's ancestor's node number in column 1. These numbers are used internally and must be unique for each node.
  The labels designate either nodes or edges. The vector \code{node.label} names internal nodes, and together with \code{tip.label}, name all nodes in the tree. The vector \code{edge.label} names all branches in the tree. All label vectors are optional, and if they are not given, internally-generated labels will be assigned. The labels, whether user-specified or internally generated, must be unique as they are used to join species data with phylogenetic trees.
  The \code{root.edge} is optional and generally not used, but provided for consistency in format (providing slots for each node in the phylogeny).   
}
\author{Ben Bolker, Thibaut Jombart}
\note{Translation functions are available from many valid tree formats. See \link{coerce-methods}.  
}
\seealso{\code{\link{coerce-methods}} for translation functions. The \linkS4class{phylo4} class, the \code{\link{check_data}}
  function to check the validity of \code{phylo4} objects. See also the \code{\link{phylo4d}} constructor, and \linkS4class{phylo4d} class.}
\examples{
# a three species tree:
mytree <- phylo4(edge=matrix(data=c(4,1, 4,5, 5,2, 5,3), ncol=2, byrow=TRUE), tip.label=c("speciesA", "speciesB", "speciesC"))
mytree
plot(mytree)

# another way to specify the same tree:
mytree <- phylo4(edge=cbind(c(4,4,5,5), c(1,5,2,3)), tip.label=c("speciesA", "speciesB", "speciesC"))

# another way:
mytree <- phylo4(edge=rbind(c(4,1), c(4,5), c(5,2), c(5,3)), tip.label=c("speciesA", "speciesB", "speciesC"))

# with branch lengths:
mytree <- phylo4(edge=rbind(c(4,1), c(4,5), c(5,2), c(5,3)), tip.label=c("speciesA", "speciesB", "speciesC"), edge.length=c(1, .2, .8, .8))
plot(mytree)
}
\keyword{classes}
