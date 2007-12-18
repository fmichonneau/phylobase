\name{subset-methods}
\docType{methods}
\alias{subset}
\alias{subset-methods}
\alias{subset,phylo-method}
\alias{subset,phylo4-method}
\alias{subset,phylo4d-method}
\title{Methods for creating subsets of phylogenies}
\description{
  Methods for creating subsets of phylogenies, based on pruning a tree to include or exclude a set of terminal taxa, to include all descendants of the MRCA of two taxa, or eventually to a subtree rooted at a given node.
}
\section{Methods}{
  \describe{
    \item{phy = "phylo4d"}{drop tips and corresponding node and tip data}
    \item{phy = "phylo4"}{drop tips}
    \item{phy = "phylo"}{drop tips}
  }
}
\note{At the moment, this simply wraps \code{ape::drop.tip}.
  Renamed from \code{drop.tip} to \code{prune} in order to
  avoid conflicts with ape.}
\usage{
\S4method{subset}{phylo4}(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,\dots)

\S4method{subset}{phylo4d}(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,\dots)

\S4method{subset}{phylo}(x,tips.include=NULL,tips.exclude=NULL,mrca=NULL,node.subtree=NULL,\dots)
}
\arguments{
  \item{x}{an object of class \code{"phylo4"}, \code{"phylo"} or \code{"phylo4d"}}
  \item{tips.include}{A vector containing tips to include in the pruned tree}
  \item{tips.exclude}{A vector containing tips to exclude from the pruned tree}
  \item{mrca}{Prune tree to include all descendents of the most recent common ancestor of the specified taxa}
  \item{node.subtree}{Not implemented yet. Evenetually will allow pruning to all descendants of a specified node}
  \item{\dots}{further arguments to be passed to other methods}
}
\details{
  To be described.

  The arguments \code{tip} can be either character or numeric. In the
  first case, it gives the labels of the tips to be deleted; in the
  second case the numbers of these labels in the vector
  \code{phy\$tip.label} are given
}
\value{
  an object of class \code{"phylo4"}
}
\author{Steven Kembel \email{skembel@berkeley.edu} and Damien de Vienne \email{damien.de-vienne@u-psud.fr}}
\examples{

}
\keyword{manip}
\keyword{methods}

