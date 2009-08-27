\name{phylo4d-hasData}

\alias{hasNodeData}
\alias{hasNodeData-methods}
\alias{hasNodeData,phylo4d-method}
\alias{hasTipData}
\alias{hasTipData-methods}
\alias{hasTipData,phylo4d-method}


\title{Tests for presence of data associated with trees stored
as phylo4d objects}
\description{
  Methods that test for the presence of data associated with trees
  stored as phylo4d objects.
}
\usage{
  \S4method{hasNodeData}{phylo4d}(x)
  \S4method{hasTipData}{phylo4d}(x)
}
\arguments{
  \item{x}{a phylo4d object}
}
\value{
  \item{\code{logical}}{return \code{TRUE} or \code{FALSE} depending
    whether data are associated with the tree (i.e., the slots
    \code{tip.data} or \code{node.data} are not empty)}
}
\details{
  The outcome of the test is based on the number of rows in the data
  frame stored in \code{tip.data} or \code{node.data}. If the number of
  rows is greater than 0, then \code{TRUE} is returned.
}

\section{Methods}{
  \describe{
    \item{hasNodeData}{\code{signature(object = "phylo4d")}: whether
      tree has internal node data}
    \item{hasTipData}{\code{signature(object = "phylo4d")}: whether tree 
      has data associated with its tips}
  }
}

\seealso{
  \code{\link{phylo4d}} constructor and \code{\linkS4class{phylo4d}} class.
}
\examples{
  data(geospiza)
  hasTipData(geospiza)  ## TRUE
  hasNodeData(geospiza) ## FALSE
}
\author{Ben Bolker, Thibault Jombart, Francois Michonneau}
\keyword{methods}
