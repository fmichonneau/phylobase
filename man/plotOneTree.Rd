% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treePlot.R
\name{plotOneTree}
\alias{plotOneTree}
\title{Plot a phylo4 object}
\usage{
plotOneTree(
  xxyy,
  type,
  show.tip.label,
  show.node.label,
  edge.color,
  node.color,
  tip.color,
  edge.width,
  rot
)
}
\arguments{
\item{xxyy}{A list created by the \code{\link{phyloXXYY}} function}

\item{type}{A character string indicating the shape of plotted tree}

\item{show.tip.label}{Logical, indicating whether tip labels should be shown}

\item{show.node.label}{Logical, indicating whether node labels should be
shown}

\item{edge.color}{A vector of colors in the order of \code{edges(phy)}}

\item{node.color}{A vector of colors indicating the colors of the node
labels}

\item{tip.color}{A vector of colors indicating the colors of the tip labels}

\item{edge.width}{A vector in the order of \code{edges(phy)} indicating the
widths of edge lines}

\item{rot}{Numeric indicating the rotation of the plot in degrees}
}
\value{
Returns no values, function invoked for the plotting side effect.
}
\description{
Plots the phylogenetic tree contained in a \code{phylo4} or \code{phylo4d}
object.
}
\examples{
library(grid)
data(geospiza)
grid.newpage()
xxyy <- phyloXXYY(geospiza)
plotOneTree(xxyy, type = 'phylogram',
  show.tip.label = TRUE, show.node.label = TRUE,
  edge.color = 'black', node.color = 'orange', tip.color = 'blue',
  edge.width = 1, rot = 0
)

grid.newpage()
pushViewport(viewport(w = 0.8, h = 0.8))
plotOneTree(xxyy, type = 'phylogram',
  show.tip.label = TRUE, show.node.label = TRUE,
  edge.color = 'black', node.color = 'orange', tip.color = 'blue',
  edge.width = 1, rot = 0
)
popViewport()

}
\seealso{
\code{treePlot}, \code{\link{phyloXXYY}}
}
\author{
Peter Cowan \email{pdc@berkeley.edu}
}
\keyword{methods}
