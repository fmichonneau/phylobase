\name{treewalk}
\alias{parent}
\alias{children}
\alias{descendants}
\alias{ancestors}
\alias{MRCA}
\alias{getnodes}
\title{tree traversal and utility functions}
\description{
  Functions for describing relationships among
  phylogenetic nodes
}
\usage{
getnodes(phy,node)
parent(phy, node)
children(phy, node)
descendants(phy, node, which=c("tips","all"))
ancestors(phy, node)
MRCA(phy,\dots)
}
\arguments{
  \item{phy}{a \code{phylo4} object (or one inheriting from
    \code{phylo4}, e.g. a \code{phylo4d} object)
  }
  \item{node}{a node number (or name)}
  \item{which}{whether to return all descendant nodes, or just tips
    (which="all" is not yet implemented)}
  \item{\dots}{a list of node numbers or names,
  or a vector of node numbers or names}
}
\value{
  \code{getNodeByLabel} takes character strings corresponding to tip or
  node labels and returns the corresponding node number;
  \code{getLabelByNode} does the opposite.
  \code{parent} and \code{children} return the numbers of the
  immediate ancestors or descendants of a node; \code{descendants}
  returns the \emph{names} of all of the descendants (tips only) of a
  node; \code{ancestors} returns the numbers of all of the ancestors of
  a node.  \code{mrca} returns the most recent common ancestor of two or
  more nodes.
}
\note{\code{MRCA} is uppercase to avoid conflict with \code{mrca}
  in ape}
\seealso{\code{mrca}, in the ape package, gives
  a list of all subtrees}
\examples{
  data(geospiza)
  getnodes(geospiza,18)
  getnodes(geospiza,"N04")
  plot(as(geospiza,"phylo4"),
    use.edge.length=FALSE,show.node.label=TRUE)
  parent(geospiza,"N11")
  children(geospiza,"N05")
  descendants(geospiza,"N11")
  ancestors(geospiza,"N11")
  MRCA(geospiza,"conirostris","difficilis","fuliginosa")
  MRCA(geospiza,"olivacea","conirostris")
}
\keyword{misc}
