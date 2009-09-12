\name{tdata}
\alias{tdata}
\alias{tdata-method}
\alias{tdata,phylo4d-method}
\alias{tdata<-}
\alias{tdata<-,phylo4d-method}
\alias{tdata<-,phylo4d,ANY-method}
\title{Retrieving or updating tip and node data in phylo4d objects}
\description{
  Method to retrieve or update tip, node or all data associated with a
  phylogenetic tree stored as a phylo4d object
}
\usage{
 \S4method{tdata}{phylo4d}(x, type=c("tip", "internal", "allnode"),
     label.type=c("row.names", "column"), empty.columns=TRUE, \dots)
 \S4method{tdata}{phylo4d,ANY}(x, type = "tip", \dots) <- value
}
\arguments{
  
    \item{x}{A \code{phylo4d} object}
    \item{type}{The type of data to retrieve or update:
      \dQuote{\code{tip}}  for data associated with tips,
      \dQuote{\code{internal}} for data  associated with internal nodes,
      \dQuote{\code{allnode}} for data associated with tip and internal
      nodes.}
    
    \item{label.type}{How the tip/node labels from the tree be
      returned? \dQuote{\code{row.names}} returns them as row names of
      the data frame, \dQuote{\code{column}} returns the labels in the
      first column of the data frame. This options is useful in the case
      of missing (\code{NA}) or non-unique labels.}
    
    \item{empty.columns}{Should columns filled with \code{NA} be
      returned?}
    
    \item{\dots}{Further arguments similar to those used by
      \code{phylo4d} (e.g. \code{match.data}), see \link{phylo4d} for
      more details.}

    \item{value}{a data frame to replace the values associated with the
      nodes specified by the argument \code{type}}

}

\section{Methods}{
  \describe{
    \item{tdata}{\code{signature(object="phylo4d")}: retrieve or update
      data associated with a tree in a \code{phylo4d} object}
  }
}
  
\value{
  \code{tdata} returns a data frame
}
\seealso{
  \code{\link{phylo4d}}
}
\examples{
   data(geospiza)
   tdata(geospiza)
   tdata(geospiza, "tip") <- 1:nTips(geospiza)
   tdata(geospiza)
   \dontshow{data(geospiza)}
}
\author{Ben Bolker, Thibaut Jombart, Francois Michonneau}
\keyword{methods}
