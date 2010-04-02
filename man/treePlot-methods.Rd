\name{treePlot-methods}
\docType{methods}
\alias{treePlot}
\alias{plot,ANY,ANY-method}
\alias{plot,pdata,missing-method}
\alias{plot,phylo4,missing-method}
\alias{treePlot-method}
\alias{treePlot,phylo4,phylo4d-method}
\title{ Phylogeny plotting }
\description{
  Plot \code{phylo4} or \code{phylo4d} objects, including associated data.
}
\section{Methods}{
  \describe{
    \item{phy = "phylo4"}{plots a tree  of class \linkS4class{phylo4}}
    \item{phy = "phylo4d"}{plots a tree with one or more
      quantitative traits contained in a \linkS4class{phylo4d} object.}
  }
}
\usage{
  \S4method{treePlot}{phylo4,phylo4d}(phy, type = c("phylogram", "cladogram", "fan"), show.tip.label = TRUE, 
  show.node.label = FALSE, tip.order = NULL, plot.data = is(phy, "phylo4d"), 
  rot = 0, tip.plot.fun = "bubbles", edge.color = "black", 
  node.color = "black", tip.color = "black", edge.width = 1, newpage = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ A \code{phylo4} or \code{phylo4d} object }
  \item{type}{ A character string indicating the shape of plotted tree }
  \item{show.tip.label}{ Logical, indicating whether tip labels should be shown }
  \item{show.node.label}{ Logical, indicating whether node labels should be shown }
  \item{tip.order}{If NULL the tree is plotted with tips in preorder, if "rev" this is reversed. Otherwise, it is a character vector of tip labels, indicating their order along the y axis (from top to bottom). Or, a numeric vector of tip node IDs indicating the order.}
  \item{plot.data}{ Logical indicating whether \code{phylo4d} data should be plotted }
  \item{rot}{ Numeric indicating the rotation of the plot in degrees }
  \item{tip.plot.fun}{ A function used to generate plot at the each tip of the phylogenetic trees }
  \item{edge.color}{ A vector of colors in the order of \code{edges(phy)} }
  \item{node.color}{ A vector of colors indicating the colors of the node labels }
  \item{tip.color}{ A vector of colors indicating the colors of the tip labels }
  \item{edge.width}{ A vector in the order of \code{edges(phy)} indicating the widths of edge lines}
  \item{newpage}{ Logical indicating whether the page should be cleared before plotting }
  \item{\dots}{ Currently unused, parameters to be passed on to \code{gpar} }  
}

\value{ No return value, function invoked for plotting side effect
}
\author{Peter Cowan \email{pdc@berkeley.edu}}

\seealso{ \code{\link{phylobubbles}}  }

\examples{

## example of plotting two grid plots on the same page
data(geospiza)
geotree <- extractTree(geospiza)
grid.newpage()
pushViewport(viewport(layout=grid.layout(nrow=1, ncol=2), name="base"))
  pushViewport(viewport(layout.pos.col=1, name="plot1"))
    treePlot(geotree, newpage=FALSE)
  popViewport()
  
  pushViewport(viewport(layout.pos.col=2, name="plot2"))
    treePlot(geotree, newpage=FALSE, rot=180)
popViewport(2)

}

\keyword{methods}
