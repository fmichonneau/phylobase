\name{phylo4-accessors}
\docType{methods}
\alias{nNodes}
\alias{nNodes-methods}
\alias{nNodes,phylo4-method}
\alias{nTips}
\alias{nTips-methods}
\alias{nTips,phylo4-method}
\alias{nTips,ANY-method}
\alias{edges}
\alias{edges-methods}
\alias{edges,phylo4-method}
\alias{nEdges}
\alias{nEdges-methods}
\alias{nEdges,phylo4-method}
\alias{edgeOrder}
\alias{edgeOrder,phylo4-method}
\alias{hasEdgeLength}
\alias{hasEdgeLength-methods}
\alias{hasEdgeLength,phylo4-method}
\alias{edgeLength}
\alias{edgeLength-methods}
\alias{edgeLength,phylo4-method}
\alias{edgeLength<-}
\alias{edgeLength<-,phylo4-method}
\alias{edgeLength<-,phylo4,ANY-method}
\alias{nodeType}
\alias{nodeType,phylo4-method}
\alias{isRooted}
\alias{isRooted-methods}
\alias{isRooted,phylo4-method}
\alias{rootEdge}
\alias{rootEdge-methods}
\alias{rootEdge,phylo4-method}
\alias{rootNode}
\alias{rootNode-methods}
\alias{rootNode,phylo4-method}
\alias{rootNode<-}
\alias{rootNode<-,phylo4-method}
\title{Methods for S4 phylogeny classes}
\description{
   Generic methods for phylogenetic trees represented as S4 classes
}
\section{Methods}{
\describe{
  \item{nTips}{\code{signature(object="phylo4")}: number of tips}
  
  \item{nNodes}{\code{signature(object="phylo4")}: number of internal nodes}
  
  \item{nEdges}{\code{signature(object = "phylo4")}: number of edges}
  
  \item{edges}{\code{signature(object = "phylo4")}: returns the edge
  matrix} 
  
  \item{hasEdgeLength}{\code{signature(object = "phylo4")}: whether
    tree has edge (branch) lengths}
  
  \item{edgeLength}{\code{signature(object = "phylo4")}: edge (branch)
    lengths (or NULL if none)}
 
  \item{nodeType}{\code{signature(object = "phylo4")}: named vector
    which has the type of node (internal, tip, root) for value, and the
    node number for name}

  \item{isRooted}{\code{signature(object = "phylo4")}: whether tree is
    rooted (i.e. has explicit root edge defined \emph{or} root node has
    <= 2 descendants)}
  
  \item{rootEdge}{\code{signature(object = "phylo4")}: root edge}
  }
}

\usage{
% what a pain: see
% https://stat.ethz.ch/pipermail/r-devel/2008-May/049653.html
\S4method{nNodes}{phylo4}(x)
\S4method{nTips}{phylo4}(x)
\S4method{edges}{phylo4}(x, order, drop.root=FALSE)
\S4method{nEdges}{phylo4}(x)
\S4method{edgeOrder}{phylo4}(x)
\S4method{hasEdgeLength}{phylo4}(x)
\S4method{edgeLength}{phylo4}(x, node)
\S4method{edgeLength}{phylo4}(x, use.names=TRUE) <- value
\S4method{nodeType}{phylo4}(phy)
\S4method{isRooted}{phylo4}(x)
\S4method{rootEdge}{phylo4}(x)
\S4method{rootNode}{phylo4}(x)
\S4method{rootNode}{phylo4}(x) <- value
}

\arguments{
  \item{x}{a phylo4/phylo4d object}
  \item{phy}{a phylo4/phylo4d object}
  \item{node}{which edge lengths to extract (indexed by descendant
    node)} 
  \item{value}{a vector of edge lengths or a node number}
  \item{order}{(not yet implemented)}
  \item{use.names}{Should the names of \code{value} be used to match
  edge lengths provided? }
  \item{drop.root}{logical: drop root row from edge matrix?}
}

\examples{
data(geospiza)
edgeLength(geospiza,5)
edgeLength(geospiza,"olivacea")
edgeLength(geospiza,5:7)
}

\keyword{methods}
