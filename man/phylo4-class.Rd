\name{phylo4-class}
\docType{class}
\alias{phylo4_orderings}
\alias{phylo-class}
\alias{phylo4-class}
\alias{$,phylo4-method}
\alias{$<-,phylo4-method}
\alias{[,phylo4,ANY,missing,missing-method}
\alias{print,phylo4-method}
\alias{show,phylo4-method}
\alias{summary,phylo4-method}
\alias{names,phylo4-method}
\title{phylo4 and extended classes}
\description{Classes for phylogenetic trees}
\section{Objects from the Class}{
  Phylogenetic tree objects can be created by calls to the \code{\link{phylo4}} constructor function.
  Translation functions from other phylogenetic packages are also available. See \code{\link{coerce-methods}}.
  }
\section{Slots}{
  \describe{
    \item{\code{edge}:}{matrix of edges}
    \item{\code{edge.length}:}{numeric vector of edge lengths}
    \item{\code{Nnode}:}{integer: number of nodes}
    \item{\code{tip.label}:}{character vector of tip labels}
    \item{\code{root.edge}:}{integer: root edge (NA if none)}
    \item{\code{order}:}{character: tree ordering (allowable values are
      listed in \code{phylo4_orderings}, currently "unknown", "preorder"
      (="pruningwise" in \code{ape}), "postorder" (="cladewise" in \code{ape}))}
  }
  }
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "phylo4")}: extract a slot}
    \item{print}{\code{signature(x = "phylo4")}: print method}
    \item{show}{\code{signature(object = "phylo4")}: show method }
    \item{summary}{\code{signature(object = "phylo4")}: summary method}
    \item{names}{\code{signature(x = "phylo4")}: gives the slot names}
  }
}
\usage{
    \S4method{print}{phylo4}(x,printlen,\dots)
    \S4method{summary}{phylo4}(x,quiet=FALSE,\dots)
    \S4method{names}{phylo4}(x)
  }
  \arguments{
    \item{x}{a phylo4 object}
    \item{printlen}{Maximum number of tip labels to print}
    \item{quiet}{a logical stating whether the results of the summary
      should be printed to the screen (FALSE, default) or not (TRUE)}
    \item{\dots}{additional arguments (for generic consistency)}
  }
  \value{
    The \code{summary} method invisibly returns a list with the
    following components: \cr
    - \code{name}: the name of the object \cr
    - \code{nb.tips}: the number of tips\cr
    - \code{nb.nodes}: the number of nodes \cr
    - \code{mean.el}: mean of edge lengths \cr
    - \code{var.el}: variance of edge lengths (estimate for population) \cr
    - \code{sumry.el}: summary (i.e. range and quartiles) of the edge
    lengths \cr
    - \code{degree}: (optional) degree (i.e. number of descendants) of
    each node; displayed only when there are polytomies \cr
    - \code{polytomy}: (optional) type of polytomy for each node:
    'node', 'terminal' (all descendants are tips) or internal (at least
    one descendant is an internal node); displayed only when there are
    polytomies \cr

   The \code{names} method returns a vector of characters corresponding
   to the names of the slots.
  }
\seealso{the \code{\link{phylo4}} constructor, the \code{\link{check_phylo4}}
  function to check the validity of \code{phylo4} objects. See also the \code{\link{phylo4d}} constructor and the \linkS4class{phylo4d} class.}
\author{Ben Bolker, Thibaut Jombart}   
\examples{
  library(ape)
  example(read.tree)
  p1 = tree.owls
  P1 = as(tree.owls,"phylo4")
  P1
  sumryP1 = summary(P1)
  sumryP1

  ## summary of a polytomous tree
  E <- matrix(c(
      8,  9, 
      9, 10, 
     10,  1, 
     10,  2, 
      9,  3, 
      9,  4, 
      8, 11, 
     11,  5, 
     11,  6, 
     11,  7, 
     NA,  8), ncol=2, byrow=TRUE) 
  
  P2 <- phylo4(E)
  plot(P2, show.node.label=TRUE)
  sumryP2 <- summary(P2)
  sumryP2
}
\keyword{classes}
