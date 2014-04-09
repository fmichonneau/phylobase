\docType{methods}
\name{treePlot-methods}
\alias{plot}
\alias{plot,ANY,ANY-method}
\alias{plot,pdata,missing-method}
\alias{plot,phylo4,missing-method}
\alias{plot,phylo4-method}
\alias{treePlot}
\alias{treePlot,phylo4,phylo4d-method}
\alias{treePlot-method}
\alias{treePlot-methods}
\title{Phylogeny plotting}
\usage{
treePlot(phy, type = c("phylogram", "cladogram", "fan"),
  show.tip.label = TRUE, show.node.label = FALSE, tip.order = NULL,
  plot.data = is(phy, "phylo4d"), rot = 0, tip.plot.fun = "bubbles",
  plot.at.tip = TRUE, edge.color = "black", node.color = "black",
  tip.color = "black", edge.width = 1, newpage = TRUE, margins = c(1.1,
  1.1, 1.1, 1.1), ...)

plot(x, y, ...)

\S4method{plot}{phylo4,missing}(x, y, ...)
}
\arguments{
  \item{phy}{A \code{phylo4} or \code{phylo4d} object}

  \item{type}{A character string indicating the shape of
  plotted tree}

  \item{show.tip.label}{Logical, indicating whether tip
  labels should be shown}

  \item{show.node.label}{Logical, indicating whether node
  labels should be shown}

  \item{tip.order}{If NULL the tree is plotted with tips in
  preorder, if "rev" this is reversed. Otherwise, it is a
  character vector of tip labels, indicating their order
  along the y axis (from top to bottom). Or, a numeric
  vector of tip node IDs indicating the order.}

  \item{plot.data}{Logical indicating whether
  \code{phylo4d} data should be plotted}

  \item{rot}{Numeric indicating the rotation of the plot in
  degrees}

  \item{tip.plot.fun}{A function used to generate plot at
  the each tip of the phylogenetic trees}

  \item{edge.color}{A vector of colors in the order of
  \code{edges(phy)}}

  \item{node.color}{A vector of colors indicating the
  colors of the node labels}

  \item{tip.color}{A vector of colors indicating the colors
  of the tip labels}

  \item{edge.width}{A vector in the order of
  \code{edges(phy)} indicating the widths of edge lines}

  \item{newpage}{Logical indicating whether the page should
  be cleared before plotting}

  \item{plot.at.tip}{should the data plots be at the tip?
  (logical)}

  \item{margins}{number of lines around the plot (similar
  to \code{par(mar)}).}

  \item{...}{Arguments to be passed to methods, such as
  \link{graphical parameters} (see \code{\link{par}}).
  Many methods will accept the following arguments:
  \describe{ \item{\code{type}}{what type of plot should be
  drawn.  Possible types are \itemize{ \item \code{"p"} for
  \bold{p}oints, \item \code{"l"} for \bold{l}ines, \item
  \code{"b"} for \bold{b}oth, \item \code{"c"} for the
  lines part alone of \code{"b"}, \item \code{"o"} for both
  \sQuote{\bold{o}verplotted}, \item \code{"h"} for
  \sQuote{\bold{h}istogram} like (or \sQuote{high-density})
  vertical lines, \item \code{"s"} for stair \bold{s}teps,
  \item \code{"S"} for other \bold{s}teps, see
  \sQuote{Details} below, \item \code{"n"} for no plotting.
  } All other \code{type}s give a warning or an error;
  using, e.g., \code{type = "punkte"} being equivalent to
  \code{type = "p"} for S compatibility.  Note that some
  methods, e.g. \code{\link{plot.factor}}, do not accept
  this.  } \item{\code{main}}{an overall title for the
  plot: see \code{\link{title}}.} \item{\code{sub}}{a sub
  title for the plot: see \code{\link{title}}.}
  \item{\code{xlab}}{a title for the x axis: see
  \code{\link{title}}.} \item{\code{ylab}}{a title for the
  y axis: see \code{\link{title}}.} \item{\code{asp}}{the
  \eqn{y/x} aspect ratio, see \code{\link{plot.window}}.} }
  }

  \item{x}{the coordinates of points in the plot.
  Alternatively, a single plotting structure, function or
  \emph{any \R object with a \code{plot} method} can be
  provided.}

  \item{y}{the y coordinates of points in the plot,
  \emph{optional} if \code{x} is an appropriate structure.}
}
\value{
No return value, function invoked for plotting side effect
}
\description{
Plot \code{phylo4} or \code{phylo4d} objects, including
associated data.
}
\section{Methods}{
  \describe{ \item{phy = "phylo4"}{plots a tree of class
  \linkS4class{phylo4}} \item{phy = "phylo4d"}{plots a tree
  with one or more quantitative traits contained in a
  \linkS4class{phylo4d} object.} }
}
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
\author{
Peter Cowan \email{pdc@berkeley.edu}
}
\seealso{
\code{\link{phylobubbles}}
}
\keyword{methods}

