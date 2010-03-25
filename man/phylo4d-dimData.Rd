\name{nData}
\alias{nData}
\alias{nData,phylo4d-method}
\title{Retrieves the number of datasets in phylo4d objects}
\description{
  Method to retrieve the number of datasets associated with a
  phylogenetic tree stored as a phylo4d object
}
\usage{

  \S4method{nData}{phylo4d}(x)

}

\arguments{
  \item{x}{A \code{phylo4d} object}
}

\details{
  \code{nData} returns the number of datasets (i.e., columns) that are
  associated with a \code{phylo4d} object.
}
\value{
  \code{nData} returns a vector.
}

\seealso{
  \code{\link{tdata}}, \code{\link{phylo4d}}
}
\examples{
  data(geospiza)
  nData(geospiza)
}
\author{
  Francois Michonnea
}
\keyword{methods}