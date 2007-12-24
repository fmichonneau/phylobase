\name{treewalk}
\alias{getAncest}
\alias{getDescend}
\alias{allDescend}
\alias{allAncest}
\alias{MRCA}
\alias{getNodeByLabel}
\title{tree traversal and utility functions}
\description{
  Functions for describing relationships among
  phylogenetic nodes
}
\usage{
getNodeByLabel(phy,lab)
getAncest(phy, node)
getDescend(phy, node)
allDescend(phy, node)
allAncest(phy, node)
MRCA(phy,\dots)
}
\arguments{
  \item{phy}{a \code{phylo4} object}
  \item{lab}{a character label}
  \item{node}{a node number (or name)}
  \item{\dots}{a list of node numbers or names,
  or a vector of node numbers or names}
}
\value{
  \code{getNodeByLabel} takes character strings corresponding to tip or
  node labels and returns the corresponding node number;
  \code{getAncest} and \code{getDescend} return the numbers of the
  immediate ancestors or descendants of a node; \code{allDescend}
  returns the \emph{names} of all of the descendants (tips only) of a
  node; \code{allAncest} returns the numbers of all of the ancestors of
  a node.  \code{mrca} returns the most recent common ancestor of two or
  more nodes.
}
\note{\code{MRCA} is uppercase to avoid conflict with \code{mrca}
  in ape}
\seealso{\code{mrca}, in the ape package, gives
  a list of all subtrees}
\examples{
  library(geiger)
  data(geospiza)
  g1 <- as(geospiza$geospiza.tree,"phylo4")
## rename nodes for clarity
  n1 <- phylobase::nTips(g1)+1
  n2 <- phylobase::nTips(g1)+nNodes(g1)
  NodeLabels(g1) <- paste("N",n1:n2,sep="")
  plot(g1,show.node.label=TRUE)
  getAncest(g1,"N19")
  getDescend(g1,"N19")
  allDescend(g1,"N19")
  allAncest(g1,"N19")
  MRCA(g1,"pallida","difficilis","fuliginosa")
  MRCA(g1,"olivacea","conirostris")
}
\keyword{misc}
