\name{tip.data.plot}
\Rdversion{1.1}
\alias{tip.data.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting trees and associated data}
\description{Plotting phylogenetic trees and associated data}
\usage{
tip.data.plot(xxyy, type = c("phylogram", "cladogram", "fan"), show.tip.label = TRUE, show.node.label = FALSE, tip.order = NULL, rot = 0, tip.plot.fun = grid.points, edge.color = "black", node.color = "black", tip.color = "black", edge.width = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xxyy}{A list created by the  \code{\link{phyloXXYY}} function}
  \item{type}{ A character string indicating the shape of plotted tree }
    \item{show.tip.label}{ Logical, indicating whether tip labels should be shown }
  \item{show.node.label}{ Logical, indicating whether node labels should be shown }
  \item{tip.order}{The order of that the tip labels should be plotted in}
  \item{rot}{ Numeric indicating the rotation of the plot in degrees }
  \item{tip.plot.fun}{A function used to plot the data elements of a \code{phylo4d} object}
  \item{edge.color}{ A vector of colors in the order of \code{edges(phy)} }
  \item{node.color}{ A vector of colors indicating the colors of the node labels }
  \item{tip.color}{ A vector of colors indicating the colors of the tip labels }
  \item{edge.width}{ A vector in the order of \code{edges(phy)} indicating the widths of edge lines}
  \item{\dots}{Additional parameters passed to tip.data.plot}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
\keyword{methods}
