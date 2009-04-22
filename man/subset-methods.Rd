\name{subset-methods}
\docType{methods}
\alias{subset}
\alias{subset-methods}
\alias{subset,phylo-method}
\alias{subset,phylo4-method}
\alias{subset,phylo4d-method}
\alias{[,phylo4-method}
\alias{[,phylo4d-method}
\title{Methods for creating subsets of phylogenies}
\description{
  Methods for creating subsets of phylogenies, based on pruning a tree
  to include or exclude a set of terminal taxa, to include all
  descendants of the MRCA of two taxa, or eventually to a subtree rooted
  at a given node.\cr

  The "[" operator can be used with \linkS4class{phylo4} and
  \linkS4class{phylo4d} objects. With a \linkS4class{phylo4} object
  'foo', \code{foo[i]} creates a subset of \code{foo} taking \code{i}
  for a tip index. With a \linkS4class{phylo4d} object
  'foo', \code{foo[i,j]} creates a subset of \code{foo} taking \code{i}
  for a tip index and \code{j} for the index of tips variables (see
  details).\cr
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
  \item{mrca}{Prune tree to include all descendants of the most recent common ancestor of the specified taxa}
  \item{node.subtree}{Prune to all descendants of a specified node}
  \item{\dots}{further arguments to be passed to other methods}
}
\details{
  The arguments \code{tip} can be either character or numeric. In the
  first case, it gives the labels of the tips to be deleted; in the
  second case the numbers of these labels in the vector
  \code{phy\$tip.label} are given.\cr

  The indices used with the \code{[} operator can be characters, integers
  or logical, just like with a vector or a matrix.\cr
  If \code{foo} is a \linkS4class{phylo4d} object, the index of
  variables ("j") is optional: \code{foo[i,j]} is the same as
  \code{foo[i]}.\cr

}
\value{
  an object of class \code{"phylo4"}
}
\author{Steven Kembel \email{skembel@berkeley.edu}\cr
  Damien de Vienne \email{damien.de-vienne@u-psud.fr}\cr
  Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}
}
\examples{
data(geospiza)

## "subset" examples
g1 <- as(geospiza,"phylo4")
nodeLabels(g1) <- paste("N", nodeId(g1), sep="")
plot(subset(g1,tips.exclude="scandens"))
plot(subset(g1,mrca=c("scandens","fortis","pauper")))
plot(subset(g1,node.subtree=18))

## "[" examples
plot(g1,show.node.label=TRUE)
plot(g1[c(1:6,14)],show.node.label=TRUE)
nodeLabels(geospiza) <- paste("N", nodeId(geospiza))
plot(geospiza[c(1:6,14)],show.node.label=TRUE)
plot(geospiza[c(1:6,14),3:5],show.node.label=TRUE)
}
\keyword{manip}
\keyword{methods}

