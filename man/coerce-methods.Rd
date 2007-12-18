\name{coerce-methods}
\docType{methods}
\alias{coerce-methods}
\alias{coerce,phylo,phylo4-method}
\alias{coerce,phylo,phylo4d-method}
\alias{coerce,multiPhylo4,multiPhylo-method}
\alias{coerce,multiPhylo,multiPhylo4-method}
\alias{coerce,phylo4,phylo-method}
\alias{coerce,phylo4d,phylo-method}
\alias{extract.tree}
\title{Converting between phylo4 and ape trees}
\description{
  Conversion functions to convert betwee \code{ape} objects
  (\code{phylo}, \code{multiPhylo}) and \code{phylo4} objects
  (\code{phylo4}, \code{phylo4d}, \code{multiPhylo4d})}
\section{Methods}{
  \describe{
    \item{from = "phylo", to = "phylo4"}{convert phylo to phylo4}
    \item{from = "phylo", to = "phylo4d"}{convert phylo to phylo4 with
      empty data}
    \item{from = "multiPhylo4", to = "multiPhylo"}{ ~~describe this method here }
    \item{from = "multiPhylo", to = "multiPhylo4"}{ ~~describe this method here }
    \item{from = "phylo4", to = "phylo"}{convert phylo4 to phylo}
    \item{from = "phylo4d", to = "phylo"}{convert phylo4d to phylo,
      dropping data (with a warning)}
  }
}
\usage{
   extract.tree(from)
}
\arguments{
  \item{from}{a \code{phylo4d} object}
  }
\author{Ben Bolker}
\examples{
    library(ape)
  example(read.tree)
## round trip conversion
  p1 = tree.owls
  P1 = as(tree.owls,"phylo4")
  P1
  identical(p1,as(P1,"phylo"))
}
\keyword{methods}
