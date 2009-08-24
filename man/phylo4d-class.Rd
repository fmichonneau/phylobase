\name{phylo4d-class}
\docType{class}
\alias{phylo4d-class}
\alias{show,phylo4d-method}
\alias{print,phylo4d-method}
\alias{head,phylo4d-method}
\alias{tail,phylo4d-method}
\alias{summary,phylo4d-method}
\alias{phylo4d,phylo4d-method}
\alias{names,phylo4d-method}
\alias{[,phylo4d,ANY,ANY,ANY-method}
\title{phylo4d class}
\description{S4 class for phylogenetic tree and data.}
\section{Objects from the Class}{
  Objects can be created from various trees and a data.frame using the
  constructor \code{phylo4d}, or using \code{new("phylo4d", \dots)} for
  empty objects.  
}
\section{Slots}{
  \describe{
    \item{\code{edge}:}{Matrix of edges}
    \item{\code{edge.label}:}{Character vector of edge (branch) labels}
    \item{\code{edge.length}:}{Numeric vector of edge (branch) lengths}
    \item{\code{Nnode}:}{Integer: number of nodes}
    \item{\code{node.data}:}{data frame for traits of internal
      nodes}
    \item{\code{node.label}:}{Character vector of node labels}
    \item{\code{tip.label}:}{Character vector of tip labels}
    % \item{\code{root.edge}:}{Integer: root edge (NA if none)}
    \item{\code{tip.data}:}{data frame for traits of tips}
    \item{\code{order}:}{character: tree ordering (allowable values are
      listed in \code{phylo4_orderings}, currently "unknown", "preorder"
      (="cladewise" in \code{ape}), and "postorder", with "cladewise"
      and "pruningwise" also allowed for compatibility with \code{ape})}
  }
}
\section{Methods}{
  \describe{    
    \item{print}{\code{signature(x = "phylo4d")}: print method}
    \item{show}{\code{signature(object = "phylo4d")}: show method }
    \item{summary}{\code{signature(object = "phylo4d")}: summary method}
    \item{names}{\code{signature(x = "phylo4d")}: gives the slots names}
    \item{head}{\code{signature(object = "phylo4d")}: show first few
nodes}
    \item{tail}{\code{signature(object = "phylo4d")}: show last few nodes}
  }
}
\usage{
    \S4method{print}{phylo4d}(x, edgeOrder = c("pretty", "real"), printall=TRUE)
    \S4method{summary}{phylo4d}(object)
    \S4method{names}{phylo4d}(x)
  }
  \arguments{
    \item{x}{a phylo4d object}
    \item{object}{a phylo4d object}
    \item{edgeOrder}{Character string indicating whether the edges should be printed
      as ordered in the tree "real" (e.g. preorder or postorder), or "pretty" printed 
      with tips collated together}
    \item{printall}{If TRUE all tip labels are printed}
  }
  \value{
    The \code{summary} method invisibly returns a list with the
    following components: \cr
    - \code{name}: the name of the object \cr
    - \code{nb.tips}: the number of tips\cr
    - \code{nb.nodes}: the number of nodes \cr
    
   The \code{names} method returns a vector of characters corresponding
   to the names of the slots.
  }
\seealso{\code{\link{coerce-methods}} for translation functions. The \code{\link{phylo4d}} constructor and the \code{\link{formatData}} function to check the validity of trees and data. See also the \code{\link{phylo4}} constructor, the \linkS4class{phylo4} class, and the \code{\link{checkPhylo4}} function to check the validity of \code{phylo4} trees. }
\author{Ben Bolker, Thibaut Jombart}   
\examples{
  library(ape)
  example(read.tree)
  obj <- phylo4d(as(tree.owls.bis,"phylo4"), data.frame(wing=1:3))
  obj
  names(obj)
  summary(obj)
}
\keyword{classes}
