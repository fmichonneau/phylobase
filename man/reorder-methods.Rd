\name{reorder-methods}
\docType{methods}
\alias{reorder-methods}
\alias{reorder,phylo-method}
\alias{reorder,phylo4-method}
\alias{reorder,phylo4d-method}
\alias{orderIndex}
\title{reordering trees within phylobase objects}
\description{
 ~~ Methods for function \code{reorder} in Package `stats' ~~
}
\section{Methods}{
  \describe{
    \item{phylo = "phylo"}{reorders a \code{phylo} object}
    \item{phylo = "phylo4"}{reorders a \linkS4class{phylo4} object}
    \item{phylo = "phylo4d"}{reorders a \linkS4class{phylo4d} object}
  }
}

\usage{
orderIndex(phy, order = 'cladewise')
\S4method{reorder}{phylo}(x, ...)
\S4method{reorder}{phylo4}(x, ...)
\S4method{reorder}{phylo4d}(x, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a \code{phylo4} or \code{phylo4d} object}
  \item{order}{(character)}
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
