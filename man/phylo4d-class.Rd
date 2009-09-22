\name{phylo4d-class}
\docType{class}
\alias{phylo4d-class}

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
    \item{\code{annote}:}{annotation data for tree (currently
      unstructured/unused by methods)}
    \item{\code{metadata}:}{metadata for node/tip data (currently
      unstructured/unused by methods)}
  }
}

% \usage{}
% \value{}
  
\seealso{
  \code{\link{coerce-methods}} for translation functions. The
  \code{\link{phylo4d}} constructor and the \code{\link{formatData}}
  function to check the validity of trees and data. See also the
  \code{\link{phylo4}} constructor, the \linkS4class{phylo4} class,
  and the \code{\link{checkPhylo4}} function to check the validity of
  \code{phylo4} trees.
}

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
