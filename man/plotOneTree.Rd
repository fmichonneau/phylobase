\name{plotOneTree}
\alias{plotOneTree}
\title{ Plot a phylo4 object }
\description{
 Plots the phylogenetic tree contained in a \code{phylo4} or \code{phylo4d} object.
}
\usage{
plotOneTree(xxyy, type, show.tip.label, show.node.label, 
  edge.color, node.color, tip.color, edge.width, rot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxyy}{ A list created by the  \code{\link{phyloXXYY}} function }
  \item{type}{ A character string indicating the shape of plotted tree }
  \item{show.tip.label}{ Logical, indicating whether tip labels should be shown }
  \item{show.node.label}{ Logical, indicating whether node labels should be shown }
  \item{edge.color}{ A vector of colors in the order of \code{edges(phy)} }
  \item{node.color}{ A vector of colors indicating the colors of the node labels }
  \item{tip.color}{ A vector of colors indicating the colors of the tip labels }
  \item{edge.width}{ A vector in the order of \code{edges(phy)} indicating the widths of edge lines}
  \item{rot}{ Numeric indicating the rotation of the plot in degrees }
}
\value{
Returns no values, function invoked for the plotting side effect.}

\author{Peter Cowan \email{pdc@berkeley.edu}}

\seealso{ \code{treePlot}, \code{\link{phyloXXYY}} }
\examples{

data(geospiza)
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

## example of plotting two grid plots on the same page
geotree <- extractTree(geospiza)
geotree2 <-
prune(geotree,descendants(geotree,MRCA(geotree,c("difficilis","scandens"))),
subtree=TRUE)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2),name="base"))
  pushViewport(viewport(layout.pos.col=1,name="plot1"))
  plot(geotree,newpage=FALSE)
  popViewport()
  pushViewport(viewport(layout.pos.col=2,name="plot2"))
  plot(geotree2,newpage=FALSE,rot=180)
  popViewport()
}

\keyword{methods}
