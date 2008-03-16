\name{pdata}
\alias{pdata}
\alias{check_pdata}
\title{Constructor for pdata (phylogenetic data) class}
\description{
  Combine data, type, comments, and metadata information
  to create a new pdata object, or check such an object
  for consistency
}
\usage{
pdata(data, type, comment, metadata)
check_pdata(object)
}
\arguments{
  \item{data}{a data frame}
  \item{type}{a factor with levels as specified by \S4class{pdata},
  the same length as \code{ncol(data)}}
  \item{comment}{a character vector, the same length as \code{ncol(data)}}
  \item{metadata}{an arbitrary list}
  \item{object}{an object of class \code{pdata}}
}
\value{
  An object of class \code{pdata}
}
\author{Ben Bolker}
\seealso{\linkS4class{pdata}}
\keyword{misc}
