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
\note{At the moment, this simply wraps \code{ape::drop.tip}.
    Renamed from \code{drop.tip} to \code{prune} in order to
    avoid conflicts with ape.}
  \usage{
\S4method{prune}{phylo4}(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, \dots)
\S4method{prune}{phylo4d}(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, \dots)
\S4method{prune}{phylo}(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0, \dots)
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
  \item{root.edge}{an integer giving the number of internal branches to
    be used to build the new root edge. This has no effect if
    \code{trim.internal = FALSE}}
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
  that have been deleted

  Note that \code{subtree = TRUE} implies \code{trim.internal = TRUE}

  To understand how the option \code{root.edge} works, see the examples
  below
}
\value{
  an object of class \code{"phylo4"}
}
\author{Emmanuel Paradis \email{Emmanuel.Paradis@mpl.ird.fr}
(original ape version)}
\examples{
require(ape)
data(bird.families)
tip <- c(
"Eopsaltriidae", "Acanthisittidae", "Pittidae", "Eurylaimidae",
"Philepittidae", "Tyrannidae", "Thamnophilidae", "Furnariidae",
"Formicariidae", "Conopophagidae", "Rhinocryptidae", "Climacteridae",
"Menuridae", "Ptilonorhynchidae", "Maluridae", "Meliphagidae",
"Pardalotidae", "Petroicidae", "Irenidae", "Orthonychidae",
"Pomatostomidae", "Laniidae", "Vireonidae", "Corvidae",
"Callaeatidae", "Picathartidae", "Bombycillidae", "Cinclidae",
"Muscicapidae", "Sturnidae", "Sittidae", "Certhiidae",
"Paridae", "Aegithalidae", "Hirundinidae", "Regulidae",
"Pycnonotidae", "Hypocoliidae", "Cisticolidae", "Zosteropidae",
"Sylviidae", "Alaudidae", "Nectariniidae", "Melanocharitidae",
"Paramythiidae","Passeridae", "Fringillidae")
plot(prune(bird.families, tip),cex=0.5)
plot(prune(bird.families, tip, trim.internal = FALSE),cex=0.5)
data(bird.orders)
plot(prune(bird.orders, 6:23, subtree = TRUE))
plot(prune(bird.orders, c(1:5, 20:23), subtree = TRUE))

### Examples of the use of `root.edge'
tr <- as(ape::read.tree(text = "(A:1,(B:1,(C:1,(D:1,E:1):1):1):1):1;"),"phylo4")
prune(tr, c("A", "B"), root.edge = 0) # = (C:1,(D:1,E:1):1);
prune(tr, c("A", "B"), root.edge = 1) # = (C:1,(D:1,E:1):1):1;
prune(tr, c("A", "B"), root.edge = 2) # = (C:1,(D:1,E:1):1):2;
prune(tr, c("A", "B"), root.edge = 3) # = (C:1,(D:1,E:1):1):3;

## Dropping tips on phylo4d objects
r1 <- rcoal(5)
d <- data.frame(a=1:5,row.names=paste("t",1:5,sep=""))
r2 <- phylo4d(r1,tip.data=d,node.data=data.frame(a=6:9))
r2 <- phylo4d(as(r1,"phylo4"),tip.data=d,node.data=data.frame(a=6:9))
nodeLabels(r2) <- as.character(sort(nodeId(r2)))
r3 <- prune(r2,1)
}
\keyword{manip}
\keyword{methods}

