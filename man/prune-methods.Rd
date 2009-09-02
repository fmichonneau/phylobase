\name{prune-methods}
\docType{methods}
\alias{prune}
\alias{prune-methods}
\alias{prune,phylo-method}
\alias{prune,phylo4-method}
\alias{prune,phylo4d-method}
\alias{na.omit,phylo4d-method}
\alias{DropTip}
\title{Methods for dropping tips}
\description{
  Methods for pruning phylogenies: remove terminal branches of
  a phylogenetic tree,
  possibly removing the corresponding internal branches.
  Preserve match to edge and node data.
}
\section{Methods}{
  \describe{
    \item{phy = "phylo4d"}{drop tips and corresponding node and tip data}
    \item{phy = "phylo4"}{drop tips}
    \item{phy = "phylo"}{drop tips}
  }
}
\note{The phylo method simply wraps \code{ape::drop.tip}.}
  \usage{
\S4method{prune}{phylo4}(phy, tip, trim.internal = TRUE, subtree = FALSE, \dots)
\S4method{prune}{phylo4d}(phy, tip, trim.internal = TRUE, subtree = FALSE, \dots)
\S4method{prune}{phylo}(phy, tip, trim.internal = TRUE, subtree = FALSE, \dots)
DropTip(phy, tip, \dots)
}
\arguments{
  \item{phy}{an object of class \code{"phylo4"}}
  \item{tip}{a vector of mode numeric or character specifying the tips
    to delete}
  \item{trim.internal}{a logical specifying whether to delete the
    corresponding internal branches}
  \item{subtree}{a logical specifying whether to output in the tree how
    many tips have been deleted and where}
  \item{\dots}{additional arguments}
}
\details{
  The argument \code{tip} can be either character or numeric. In the
  first case, it gives the labels of the tips to be deleted; in the
  second case the numbers of these labels in the vector
  \code{tipLabels(phy)} are given

  If \code{trim.internal = FALSE}, the new tips are given \code{NA} as
  labels, unless there are node labels in the tree in which case they
  are used

  If \code{subtree = TRUE}, the returned tree has one or several
  terminal branches indicating how many tips have been removed (with a
  label \code{"[x_tips]"}). This is done for as many monophyletic groups
  that have been deleted. This option is only implemented for the
  phylo method, and thus if used on phylo4 or phylo4d objects, they will
  be coerced to phylo

  Note that \code{subtree = TRUE} implies \code{trim.internal = TRUE}

  Handling of root edge lengths differs between \code{phylobase} and
  \code{ape}. For phylo4/phylo4d objects, the new root edge length
  will be the total distance back through the original root edge (or NA
  if the original edge had length NA)

  A \code{root.edge} argument can be passed to the phylo method. To
  understand how this option works, see the examples below
}
\value{
  an object of class \code{"phylo4"}
}
\author{
  Jim Regetz \email{regetz@nceas.ucsb.edu} (phylo4 method),
  Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr
    (original ape version used in phylo method)}
}
\examples{
data(geospiza)

## Pruning phylo4 objects
geotree <- extractTree(geospiza)
tip <- c("difficilis", "fortis", "fuliginosa", "fusca", "olivacea",
    "pallida", "parvulus", "scandens")
plot(prune(geotree, tip))
plot(prune(geotree, tip, trim.internal = FALSE))
plot(prune(geotree, 1:11, subtree = TRUE))

## Pruning phylo4d objects
r1 <- as(rcoal(5), "phylo4")
td <- data.frame(a=1:5, row.names=paste("t", 1:5, sep=""))
nd <- data.frame(a=6:9, row.names=nodeId(r1, "internal"))
r2 <- phylo4d(r1, tip.data=td, node.data=nd)
r3 <- prune(r2, "t1")

## Pruning phylo objects, including the use of `root.edge'
require(ape)
tr <- as(ape::read.tree(text = "(A:1,(B:1,(C:1,(D:1,E:1):1):1):1):1;"),"phylo4")
prune(tr, c("A", "B"), root.edge = 0) # = (C:1,(D:1,E:1):1);
prune(tr, c("A", "B"), root.edge = 1) # = (C:1,(D:1,E:1):1):1;
prune(tr, c("A", "B"), root.edge = 2) # = (C:1,(D:1,E:1):1):2;
prune(tr, c("A", "B"), root.edge = 3) # = (C:1,(D:1,E:1):1):3;

}
\keyword{manip}
\keyword{methods}

