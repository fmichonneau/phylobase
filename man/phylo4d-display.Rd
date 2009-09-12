\name{phylo4d-display}
\docType{methods}
\alias{show,phylo4d-method}
\alias{print,phylo4d-method}
\alias{head,phylo4d-method}
\alias{tail,phylo4d-method}
\alias{summary,phylo4d-method}
\alias{names,phylo4d-method}

\title{Display methods for phylo4d objects}

\description{
  Methods used to display information about the data and the tree for
  phylo4d objects.
}

\section{Methods}{
  \describe{    
    \item{print}{\code{signature(x = "phylo4d")}: print method}
    \item{show}{\code{signature(object = "phylo4d")}: show method }
    \item{summary}{\code{signature(object = "phylo4d")}: summary method}
    \item{names}{\code{signature(x = "phylo4d")}: gives the slots names}
    \item{head}{\code{signature(object = "phylo4d")}: show first few
      nodes}
    \item{tail}{\code{signature(object = "phylo4d")}: show last few
      nodes}
  }
}

\usage{
    \S4method{print}{phylo4d}(x, edgeOrder = c("pretty", "real"),
    printall=TRUE) 
    \S4method{summary}{phylo4d}(object, quiet=FALSE)
    \S4method{names}{phylo4d}(x)
}

\arguments{
  \item{x}{a phylo4d object}
  \item{object}{a phylo4d object}
  \item{edgeOrder}{Character string indicating whether the edges should be printed
    as ordered in the tree "real" (e.g. preorder or postorder), or "pretty" printed 
    with tips collated together}
  \item{quiet}{Should the summary be displayed on screen?}
  \item{printall}{If TRUE all tip labels are printed}
}
\value{
    The \code{summary} method invisibly returns a list with the
    following components:
    \item{\code{name}}{the name of the object}
    \item{\code{nb.tips}}{the number of tips}
    \item{\code{nb.nodes}}{the number of nodes}
    \item{\code{mean.el}}{mean of edge lengths}
    \item{\code{var.el}}{variance of edge lengths (estimate for
      population) }
    \item{\code{sumry.el}}{summary (i.e. range and quartiles) of the
      edge lengths}
    \item{\code{degree}}{(optional) type of polytomy for each node:
      \sQuote{node}, \sQuote{terminal} (all descendants are tips) or 
      \sQuote{internal} (at least one descendant is an internal node);
      displayed only when there are polytomies}
    \item{\code{sumry.tips}}{(optional) summary for the data associated
      with the tips}
    \item{\code{sumry.nodes}}{(optional) summary for the data associated
      with the internal nodes}    
    
    The \code{names} method returns a vector of characters corresponding
    to the names of the slots.
}

\seealso{
  \code{\link{phylo4d}} constructor and \code{\linkS4class{phylo4d}} class.
}

\author{Ben Bolker, Thibaut Jombart}
% \examples{}
\keyword{methods}
  