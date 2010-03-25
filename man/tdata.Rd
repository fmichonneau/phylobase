\name{tdata}
\alias{tdata}
\alias{tdata-method}
\alias{tdata,phylo4d-method}
\alias{tdata<-}
\alias{tdata<-,phylo4d-method}
\alias{tdata<-,phylo4d,ANY-method}
\alias{tipData}
\alias{tipData-method}
\alias{tipData,phylo4d-method}
\alias{tipData<-}
\alias{tipData<-,phylo4d-method}
\alias{tipData<-,phylo4d,ANY-method}
\alias{nodeData}
\alias{nodeData-method}
\alias{nodeData,phylo4d-method}
\alias{nodeData<-}
\alias{nodeData<-,phylo4d-method}
\alias{nodeData<-,phylo4d,ANY-method}
\title{Retrieving or updating tip and node data in phylo4d objects}
\description{
  Methods to retrieve or update tip, node or all data associated with a
  phylogenetic tree stored as a phylo4d object
}
\usage{
  \S4method{tdata}{phylo4d}(x, type=c("all", "tip", "internal"),
    label.type=c("row.names", "column"), empty.columns=TRUE)
  \S4method{tdata}{phylo4d,ANY}(x, type=c("all", "tip", "internal"),
    merge.data=TRUE, clear.all=FALSE, \dots) <- value
  \S4method{tipData}{phylo4d}(x, \dots)
  \S4method{tipData}{phylo4d,ANY}(x, \dots) <- value
  \S4method{nodeData}{phylo4d}(x, \dots)
  \S4method{nodeData}{phylo4d,ANY}(x, \dots) <- value
}
\arguments{
  
    \item{x}{A \code{phylo4d} object}

    \item{type}{The type of data to retrieve or update:
      \dQuote{\code{all}} (default) for data associated with both tip
      and internal nodes, \dQuote{\code{tip}} for data associated with
      tips only, \dQuote{\code{internal}} for data associated with
      internal nodes only.}
    
    \item{label.type}{How should the tip/node labels from the tree be
      returned? \dQuote{\code{row.names}} returns them as row names of
      the data frame, \dQuote{\code{column}} returns them in the first
      column of the data frame. This options is useful in the case of
      missing (\code{NA}) or non-unique labels.}
    
    \item{empty.columns}{Should columns filled with \code{NA} be
      returned?}

    \item{merge.data}{if tip or internal node data are provided and data
      already exists for the other type, this determines whether columns
      with common names will be merged together (default TRUE). If
      FALSE, columns with common names will be preserved separately,
      with \dQuote{.tip} and \dQuote{.node} appended to the names. This
      argument has no effect if tip and node data have no column names
      in common, or if type=\dQuote{all}.}
    
    \item{clear.all}{If only tip or internal node data are to be
      replaced, should data of the other type be dropped?}

    \item{\dots}{For the tipData and nodeData accessors, further
      arguments to be used by tdata. For the replacement forms, further
      arguments to be used by \code{formatData} (e.g.
      \code{match.data}), see \link{formatData} for more details.}

    \item{value}{a data frame (or object to be coerced to one) to
      replace the values associated with the nodes specified by the
      argument \code{type}}

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
   tipData(geospiza) <- 1:nTips(geospiza)
   tdata(geospiza)
   \dontshow{data(geospiza)}
}
\author{Ben Bolker, Thibaut Jombart, Francois Michonneau}
\keyword{methods}
