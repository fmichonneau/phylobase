\name{phylo4-class}
\docType{class}
\alias{phylo4}
\alias{phylo-class}
\alias{phylo4-class}
\alias{multiPhylo-class}
\alias{multiPhylo4-class}
\alias{$,phylo4-method}
\alias{print,phylo4-method}
\alias{show,phylo4-method}
\alias{summary,phylo4-method}
\alias{names,phylo4-method}
\title{phylo4 and extended classes}
\description{Classes for phylogenetic trees}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("phylo4", ...)}
  Translation from other phylogenetic packages ...
  }
\section{Slots}{
  \describe{
    \item{\code{edge}:}{Matrix of edges}
    \item{\code{edge.length}:}{Numeric vector of edge lengths}
    \item{\code{Nnode}:}{Integer: number of nodes}
    \item{\code{tip.label}:}{Character vector of tip labels}
    \item{\code{root.edge}:}{Integer: root edge (NA if none)}
    \item{\code{phylolist}:}{(multiPhylo4) list of phylogenetic trees (each of which
    may be of type \code{phylo4d}, containing node and/or edge data)}
}
}
\section{Methods}{
  \describe{
    \item{\$}{\code{signature(x = "phylo4")}: extract a slot}
    \item{print}{\code{signature(x = "phylo4")}: print method}
    \item{show}{\code{signature(object = "phylo4")}: show method }
    \item{summary}{\code{signature(object = "phylo4")}: summary method}
    \item{names}{\code{signature(x = "phylo4")}: gives the slots names}
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

   The \code{names} method returns a vector of characters corresponding
   to the names of the slots.
  }
  \author{Ben Bolker, Thibaut Jombart}   
\examples{
  library(ape)
  example(read.tree)
  p1 = tree.owls
  P1 = as(tree.owls,"phylo4")
  P1
  sumryP1 = summary(P1)
  sumryP1
}
\keyword{classes}
