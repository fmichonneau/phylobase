\name{plot for phylo4 and phylo4d}
\docType{methods}
\alias{plot,phylo4,missing-method}
\alias{plot,phylo4d,missing-method}
\title{Plot a phylogenetic tree alone or with data}
\description{
  These methods display graphics for objects of class \code{phylo4} and
  \linkS4class{phylo4d}.
  Plotting of phylogenies calls \code{plot.phylo} from the \code{ape}
  package.
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
\S4method{plot}{phylo4d,missing}(x, treetype=c("phylogram","cladogram"), symbol=c("circles", "squares"), center=TRUE, scale=TRUE, legend=TRUE, grid=TRUE, box=TRUE, show.tip.label=TRUE, show.node.label=TRUE, show.var.label=TRUE, ratio.tree=1/3, font=3, tip.label=x@tip.label, var.label=colnames(x@tip.data), cex.symbol=1, cex.label=1, cex.legend=1, \dots)

}
\arguments{
  \item{x}{an object of class \linkS4class{phylo4} or
    \linkS4class{phylo4d}}
  \item{treetype}{the type of tree to be plotted ("phylogram" or "cladogram")}
  \item{symbol}{the type of symbol used to represent data ("circles" or "squares")}
  \item{center}{a logical stating whether variables should be centred
    (TRUE, default) or not (FALSE)}
  \item{scale}{a logical stating whether variables should be scaled
    (TRUE, default) or not (FALSE)}
  \item{legend}{a logical stating whether a legend should be added to
    the plot (TRUE) or not (FALSE, default)}
  \item{grid}{a logical stating whether a grid should be added to the plot
    (TRUE, default) or not (FALSE)}
  \item{box}{a logical stating whether a box should be added around the plot
    (TRUE, default) or not (FALSE)}
 \item{show.tip.label}{a logical stating whether tip labels should be
   printed (TRUE, default) or not (FALSE)}
 \item{show.node.label}{a logical stating whether node labels should be
   printed (TRUE, default) or not (FALSE)}
 \item{show.var.label}{a logical stating whether labels of variables should be
   printed (TRUE, default) or not (FALSE)}
 \item{ratio.tree}{the proportion of width of the figure occupied by the tree}
 \item{font}{an integer specifying the type of font for the labels: 1
   (plain text), 2 (bold), 3 (italic, default), or 4 (bold italic).}
 \item{tip.label}{a character vector giving the tip labels}
 \item{var.label}{a character vector giving the labels of variables}
 \item{cex.symbol}{a numeric giving the factor scaling the symbols}
 \item{cex.label}{a numeric giving the factor scaling all labels}
 \item{cex.legend}{a numeric giving the factor scaling the legend}
 \item{\dots}{further arguments to be passed to plot methods from \code{ape}}
}
\seealso{the  \linkS4class{phylo4} and \linkS4class{phylo4d} classes;
  \code{\link[pkg:ape]{plot.phylo}} from the \code{ape} package,
    \code{\link[pkg:ade4]{symbols.phylog}},
    \code{\link[pkg:ade4]{table.phylog}} and
    \code{\link[pkg:ade4]{dotchart.phylog}}.}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\examples{
if(require(ape)){
## build data
example("phylo4d")
obj1 <- obj2 <- obj3 <- phylo4d(as(tree.owls,"phylo4"),data.frame(wing=1:4,color=factor(c("b","w","b","b")), tail=runif(4)*10), use.tip.names=FALSE)

obj2@tip.data <- as.data.frame(obj2@tip.data[,1])
obj3@tip.data <- cbind(obj1@tip.data,obj2@tip.data)
obj4 <- obj1
obj4$tip.data[2,3] <- NA
obj4$tip.data[1,1] <- NA

## several plots
plot(obj1)
plot(obj2,box=FALSE,symbol="sq",var.lab="My legend for the variable",ratio=.6,center=FALSE,scale=FALSE,leg=TRUE,cex.leg=1.5)
plot(obj3,leg=TRUE,cex.leg=1.2,treetype="clado",adj=-.5,srt=90)
plot(obj4,edge.color=rainbow(6),edge.width=1:6)

if(require(ade4)){
data(mjrochet)
temp <- as(read.tree(text=mjrochet$tre),"phylo4")
obj <- phylo4d(x=temp,tip.data=mjrochet$tab)
obj
plot(obj,cex.lab=.5,show.node=FALSE,cex.sym=.6)
}
}
}
\keyword{misc}
\keyword{hplot}