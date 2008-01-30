\name{check_phylo4}
\alias{check_phylo4}
\title{Validity checking for phylo4 objects}
\description{
  Basic checks on the validity of S4 phylogenetic objects
}
\usage{
check_phylo4(object)
}
\arguments{
  \item{object}{A prospective S4 object}
}
\value{
  As required by \code{\link{validObject}}, returns an
  error string (describing problems) or TRUE if everything is OK
}
\seealso{
  the \code{\link{phylo4}} constructor, the \linkS4class{phylo4} class. See also the \code{\link{check_data}}, the \code{\link{phylo4d}} constructor and the \linkS4class{phylo4d} class. See \code{\link{coerce-methods}} for translation functions.
}
\author{Ben Bolker}
\keyword{misc}
