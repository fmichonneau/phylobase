\name{phylobubbles}
\alias{phylobubbles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bubble plots for phylo4d objects }
\description{
Plots either circles or squares corresponding to the magnitude of each cell of a \code{phylo4d} object.}
\usage{
phylobubbles(type, place.tip.label, show.node.label, tip.order, rot, edge.color, node.color, tip.color, edge.width, newpage, ..., XXYY, square = FALSE, grid = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{the type of plot }
  \item{place.tip.label}{A string indicating whether labels should be plotted to the right or to the left of the bubble plot}
  \item{show.node.label}{A logical indicating whether internal node labels should be plotted }
  \item{tip.order}{The order of that the tip labels should be plotted in}
  \item{rot}{The number of degrees that the plot should be rotated }
  \item{edge.color}{A vector of colors for the tree edge segments}
  \item{node.color}{A vector of colors for the coloring the nodes}
  \item{tip.color}{A vector of colors for the coloring the tip labels}
  \item{edge.width}{A vector of line widths for the tree edges}
  \item{newpage}{Logical to control whether the device is cleared before plotting, useful for adding plot inside other plots}
  \item{\dots}{Additional parameters passed to the bubble plotting functions}
  \item{XXYY}{The out put from the phyloXXYY function}
  \item{square}{Logical indicating whether the plot 'bubbles' should be squares}
  \item{grid}{A logical indicating whether a grey grid should be plotted behind the bubbles}

}

\author{Peter Cowan \email{pdc@berkeley.edu}}

\seealso{\code{\link{phyloXXYY}}, \code{treePlot} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{methods}