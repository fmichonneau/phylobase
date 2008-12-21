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
    \item{\code{edge.length}:}{Numeric vector of edge lengths}
    \item{\code{Nnode}:}{Integer: number of nodes}
    \item{\code{tip.label}:}{Character vector of tip labels}
    \item{\code{root.edge}:}{Integer: root edge (NA if none)}
    \item{\code{tipdata}:}{data frame for traits of tips}
    \item{\code{nodedata}:}{data frame for traits of internal
      nodes}
    \item{\code{edgedata}:}{data frame for traits of edges}
  }
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "phylo4d")}: extract a slot}
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
    \S4method{print}{phylo4d}(x,printlen,\dots)
    \S4method{summary}{phylo4d}(x,\dots)
    \S4method{names}{phylo4d}(x)
  }
  \arguments{
    \item{x}{a phylo4d object}
    \item{printlen}{maximum number of tip labels to print}
    \item{\dots}{additional arguments (for generic consistency)}
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
\seealso{\code{\link{coerce-methods}} for translation functions. The \code{\link{phylo4d}} constructor and the \code{\link{check_data}} function to check the validity of trees and data. See also the \code{\link{phylo4}} constructor, the \linkS4class{phylo4} class, and the \code{\link{check_phylo4}} function to check the validity of \code{phylo4} trees. }
\author{Ben Bolker, Thibaut Jombart}   
\examples{
  library(ape)
  example(read.tree)
  obj <- phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3), use.tip.names=FALSE)
  obj
  names(obj)
  summary(obj)
}
\keyword{classes}
