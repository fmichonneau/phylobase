\name{plot for phylo4 and phylo4d}
\docType{methods}
\alias{plot,phylo4,missing-method}
\alias{plot,phylo4d,missing-method}
\title{Plot a phylogenetic tree alone or with data}
\description{
  These methods display graphics for objects of class \code{phylo4} and
  \linkS4class{phylo4d}.
  They are currently wrappers to \code{plot.phylo} from the \code{ape}
  package (for \linkS4class{phylo4} objects) and to various graphics from the
  \code{ade4} package (for \linkS4class{phylo4d} objects).
}
\section{Methods}{
  \describe{
\item{x = "phylo4", y = "missing"}{plots a tree  of class \linkS4class{phylo4}}

\item{x = "phylo4d", y = "missing"}{plots a tree with one or more
  quantitative traits contained in a \linkS4class{phylo4d} object.}

   }
 }
 \usage{
\S4method{plot}{phylo4,missing}(x, \dots)
\S4method{plot}{phylo4d,missing}(x, type=c("symbols", "squares",
"dotchart"), show.tip.label=TRUE, show.node.label=FALSE,
cex.label=par("cex"), cex.symbol=1, \dots)

}
\arguments{
  \item{x}{an object of class \linkS4class{phylo4} or \linkS4class{phylo4d}}
  \item{type}{a character string giving the type of graphic to be drawn:
  "symbols" is used for a single variable; "squares" and "dotchart" are
  used for multivariate data.}
\item{show.tip.label}{a logical stating whether tip labels should be
  printed (TRUE, default) or not (FALSE)}
\item{show.node.label}{a logical stating whether node labels should be
  printed (TRUE) or not (FALSE, default)}
\item{cex.label}{a numeric giving the factor scaling the tip and node labels}
\item{cex.symbol}{a numeric giving the factor scaling the symbols (not
  used if type=="dotchart")}
  \item{\dots}{further arguments to be passed to plot methods from \code{ape}
    or \code{ade4}.}
}
\seealso{the  \linkS4class{phylo4} and \linkS4class{phylo4d} classes;
  \code{\link[pkg:ape]{plot.phylo}} from the \code{ape} package,
    \code{\link[pkg:ade4]{symbols.phylog}},
    \code{\link[pkg:ade4]{table.phylog}} and \code{\link[pkg:ade4]{dotchart.phylog}}.}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\examples{
if(require(ape) & require(ade4)){
example(read.tree)
obj1 <- phylo4d(as(tree.owls,"phylo4"),data.frame(wing=1:4), use.tip.name=FALSE)
obj2 <- phylo4d(as(tree.owls,"phylo4"),data.frame(wing=1:4, tail=4:1, body.length=c(10,15,5,20)), use.tip.name=FALSE)

par(mfrow=c(2,2))
plot(as(obj1,"phylo4"))
plot(obj1)
plot(obj2, show.node=TRUE)
plot(obj2, type="dotchart")

}
}
\keyword{misc}
