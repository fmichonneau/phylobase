\name{phyloXXYY}
\alias{phyloXXYY}
\title{ Calculate node x and y coordinates }
\description{
Calculates the node x and y locations for plotting a phylogenetic tree.}
\usage{
phyloXXYY(phy, tip.order = NULL)
}
\arguments{
  \item{phy}{ A \code{phylo4} or \code{phylo4d} object. }
  \item{tip.order}{ Unused. A character vector indicating the ordering of tip nodes}
}
\details{
The y coordinates of the tips are evenly spaced from 0 to 1 in pruningwise order.  Ancestor y nodes are given the mean value of immediate descendants.  The root is given the x coordinate 0 and descendant nodes are placed according to the cumulative branch length from the root, with a maximum x value of 1.
}
\value{
  \item{yy}{Internal node and tip y coordinates}
  \item{xx}{Internal node and tip x coordinates}
  \item{traverse}{The tree traversal order}
  \item{phy}{A \code{phylo4} or \code{phylo4d} object}
  \item{phy.orig}{The original \code{phylo4} or \code{phylo4d} object}
}

\author{Peter Cowan \email{pdc@berkeley.edu}}

\seealso{\code{treePlot}, \code{\link{plotOneTree}}}
\examples{

data(geospiza)
coor <- phyloXXYY(geospiza)
plot(coor$xx, coor$yy, pch = 20)

}
\keyword{methods}
