\name{checkPhylo4}
\alias{checkPhylo4}
\alias{checkTree}
\title{Validity checking for phylo4 objects}
\description{
  Basic checks on the validity of S4 phylogenetic objects
}
\usage{
checkPhylo4(object)
checkTree(object,warn="retic",err=NULL)
}
\arguments{
  \item{object}{A prospective S4 object}
  \item{warn}{a character vector listing phenomena to warn about:
    current options are "poly" and "retic"}
  \item{err}{a character vector listing phenomena to trigger errors}
}
\value{
  As required by \code{\link{validObject}}, returns an
  error string (describing problems) or TRUE if everything is OK
}
\note{\code{checkPhylo4} is an (inflexible) wrapper for
  \code{checkTree}, which uses the default settings.
  The rules for \code{phylo4} objects essentially follow
  those for \code{phylo} objects from the \code{ape} package,
  which are in turn defined in
  http://ape.mpl.ird.fr/misc/FormatTreeR_4Dec2006.pdf.
  These are essentially that:
  \itemize{
    \item if the tree has edge lengths defined, the number of edge lengths must match the number of edges;
    \item the number of tip labels must match the number of tips;
    \item in a tree with \code{ntips} tips and \code{nnodes} (total) nodes, nodes 1 to \code{ntips} must be tips
    \item if the tree is rooted, the root must be node number
    \code{ntips+1} and the root node must be the first row of the edge matrix
  }
}
\seealso{
  the \code{\link{phylo4}} constructor and \linkS4class{phylo4} class;
  \code{\link{checkData}}, the \code{\link{phylo4d}} constructor and
  the \linkS4class{phylo4d}
  class do checks for the data associated with trees.
  See \code{\link{coerce-methods}} for translation functions.
}
\author{Ben Bolker}
\keyword{misc}
