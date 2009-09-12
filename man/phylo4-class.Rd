\name{phylo4-class}
\docType{class}
\alias{phylo4_orderings}
\alias{phylo-class}
\alias{phylo4-class}
\title{The phylo4 class}
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
      listed in \code{phylo4_orderings}, currently "unknown",
      "preorder" (="cladewise" in \code{ape}), and "postorder",
      with "pruningwise" and "cladewise" also allowed for compatibility with
      \code{ape})}
      \item{\code{annote}:}{annotation data for tree (currently
      unstructured/unused by methods)}
    }
}

% \usage{}
% \arguments{}
% \value{}
\seealso{
  The \code{\link{phylo4}} constructor, the \code{\link{checkPhylo4}}
  function to check the validity of \code{phylo4} objects. See also the
  \code{\link{phylo4d}} constructor and the \linkS4class{phylo4d}
  class.} 

\author{Ben Bolker, Thibaut Jombart}
% \examples{}
\keyword{classes}