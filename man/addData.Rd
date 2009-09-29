\name{addData}
\alias{addData}
\alias{addData-methods}
\alias{addData,phylo4-method}
\alias{addData,phylo4d-method}

\title{Adding data to a phylo4 or a phylo4d object}
\description{
  \code{addData} adds data to a \code{phylo4} (converting it in a
  \code{phylo4d} object) or to a \code{phylo4d} object
}
\usage{
\S4method{addData}{phylo4}(x, tip.data, node.data, all.data,
pos=c("after", "before"), merge.data=TRUE, match.data=TRUE, \dots)
\S4method{addData}{phylo4d}(x, tip.data, node.data, all.data,
pos=c("after", "before"), merge.data=TRUE, match.data=TRUE, \dots)
}
\arguments{
  \item{x}{a phylo4 or a phylo4d object}
  \item{tip.data}{a data frame for tip data}
  \item{node.data}{a data frame for node data}
  \item{all.data}{a data frame for both tip an node data}
  \item{match.data}{(logical), should the row names provided in
    \code{tip.data}, \code{node.data} and/or \code{all.data} be used to
    match data against the tree labels?}
  \item{pos}{should the new data provided be bound \code{before} or
    \code{after} the pre-existing data?}
  \item{merge.data}{if both \code{tip.data} and \code{node.data} are
    provided, it determines whether columns with common names will be
    merged together (default TRUE). If FALSE, columns with common names
    will be preserved separately, with ".tip" and ".node" appended to
    the names. This argument has no effect if \code{tip.data} and
    \code{node.data} have no column names in common.}
  \item{\dots}{additional arguments to be passed to \link{formatData}}
}
\value{
  \code{addData} returns a \code{phylo4d} object. 
}
\details{
  If any column names in the original data are the same as columns in
  the new data, ".old" is appended to the former column names and ".new"
  is appended to the new column names.

  The option \code{pos} is ignored (silently) if \code{x} is a
  \code{phylo4} object. It is provided for compatibility reasons.
}
\seealso{
  \code{\link{tdata}} for extracting or updating data and
  \code{\link{phylo4d}} constructor.
}
\examples{
  data(geospiza)
  nDt <- data.frame(a=rnorm(nNodes(geospiza)), b=1:nNodes(geospiza),
    row.names=nodeId(geospiza, "internal"))
  t1 <-addData(geospiza, node.data=nDt)
}
\author{Francois Michonneau}
\keyword{methods}
