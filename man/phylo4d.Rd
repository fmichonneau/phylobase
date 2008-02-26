\name{phylo4d-methods}
\docType{methods}
\alias{phylo4d}
\alias{phylo4d-methods}
\alias{phylo4d,phylo4-method}
\alias{phylo4d,matrix-method}
\alias{phylo4d,phylo-method}
\title{Combine a phylogenetic tree with data}
\description{
  \code{phylo4d} is a generic constructor which merges a phylogenetic tree with
  a data frame to create a combined object of class \code{phylo4d}
}
\section{Methods}{
  \describe{
\item{x = "phylo4"}{merges a tree of class \code{phylo4} with a
  data.frame into a \code{phylo4d} object}

\item{x = "matrix"}{merges a matrix of tree edges similar to the edge
  slot of a \code{phylo4} object (or to \$edge of a \code{phylo} object) with a
  data.frame into a \code{phylo4d} object}

\item{x = "phylo"}{merges a tree of class \code{phylo} with a
  data.frame into a \code{phylo4d} object }    
   }
 }
 \usage{
\S4method{phylo4d}{phylo4}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL,...)
\S4method{phylo4d}{phylo}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL,...)
\S4method{phylo4d}{matrix}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL,...)
}
\arguments{
  \item{x}{an object of class \code{phylo4}, \code{phylo} or a matrix of
  edges (see above)}
\item{tip.data}{a data frame for tips data}
\item{node.data}{a data frame for nodes data}
\item{all.data}{a data frame for all (i.e. tips and nodes) data. In such
  case, first rows should correspond to tips, last rows to nodes.}
\item{\dots}{further arguments to be passed to other methods. Used to
  control the validity checking of the created object (see \code{\link{check_data}}).}
}
\value{
  An object of class \linkS4class{phylo4d}.
}
\seealso{
\code{\link{coerce-methods}} for translation functions. The \linkS4class{phylo4d} class, the \code{\link{check_data}}
  function to check the validity of \code{phylo4d} objects; \linkS4class{phylo4} class and \link{phylo4} constructor.}
\author{Ben Bolker, Thibaut Jombart}
\note{
  Checking on matches will be done by the validity checker
  (label matches between data and tree tips, number of rows
  of data vs. number of nodes/tips/etc.)
}
\examples{
library(ape)
example(read.tree)
try(phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3)), silent=TRUE)
obj <- phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3),
use.tip.names=FALSE)
obj
print(obj)
}
\keyword{misc}
