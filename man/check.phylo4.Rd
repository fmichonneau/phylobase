\name{check_phylo4}
\alias{check_phylo4}
\alias{check_tree}
\title{Validity checking for phylo4 objects}
\description{
  Basic checks on the validity of S4 phylogenetic objects
}
\usage{
check_phylo4(object)
check_tree(object,warn="retic",err=NULL)
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
\note{\code{check_phylo4} is an (inflexible) wrapper for
  \code{check_tree} that uses the default settings}
\seealso{
  the \code{\link{phylo4}} constructor and \linkS4class{phylo4} class;
  \code{\link{check_data}}, the \code{\link{phylo4d}} constructor and
  the \linkS4class{phylo4d}
  class do checks for the data associated with trees.
  See \code{\link{coerce-methods}} for translation functions.
}
\author{Ben Bolker}
\keyword{misc}
