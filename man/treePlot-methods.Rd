\name{treePlot-methods}
\docType{methods}
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
  rot = 0, tip.plot.fun = "bubbles", plot_by_tip = TRUE, edge.color = "black", 
  node.color = "black", tip.color = "black", edge.width = 1, newpage = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{ A \code{phylo4} or \code{phylo4d} object }
  \item{type}{ A character string indicating the shape of plotted tree }
  \item{show.tip.label}{ Logical, indicating whether tip labels should be shown }
  \item{show.node.label}{ Logical, indicating whether node labels should be shown }
  \item{tip.order}{ Unused, a character vector of tip labels, indicating they're order along the y axis }
  \item{plot.data}{ Logical indicating whether \code{phylo4d} data should be plotted }
  \item{rot}{ Numeric indicating the rotation of the plot in degrees }
  \item{tip.plot.fun}{ A function used to generate plot at the each tip of the phylogenetic trees }
  \item{plot_by_tip} {TODO} % TODO what does the plot_by_tip argument do
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


\keyword{methods}
