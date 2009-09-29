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
    \item{\code{edge}:}{Matrix of edges}
    \item{\code{edge.label}:}{Character vector of edge (branch) labels}
    \item{\code{edge.length}:}{Numeric vector of edge (branch) lengths}
    \item{\code{label}:}{Character vector of tip (and optionally
      internal) node labels}
    \item{\code{order}:}{character: tree ordering (allowable values are
      listed in \code{phylo4_orderings}, currently "unknown", "preorder"
      (="cladewise" in \code{ape}), and "postorder", with "cladewise"
      and "pruningwise" also allowed for compatibility with \code{ape})}
    \item{\code{annote}:}{annotation data for tree (currently
      unstructured/unused by methods)}
    \item{\code{metadata}:}{metadata for node/tip data (currently
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
