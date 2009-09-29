\name{phylo4-display}
\docType{method}
\alias{print,phylo4-method}
\alias{show,phylo4-method}
\alias{head,phylo4-method}
\alias{tail,phylo4-method}
\alias{summary,phylo4-method}
\alias{names,phylo4-method}
\title{Displaying phylo4 object}
\description{Display methods for phylo4 and phylo4d phylogenetic trees}


\section{Methods}{
  \describe{    
    \item{print}{\code{signature(x = "phylo4")}: print method}
    \item{show}{\code{signature(object = "phylo4")}: show method }
    \item{summary}{\code{signature(object = "phylo4")}: summary method}
    \item{names}{\code{signature(x = "phylo4")}: gives the slot names}
    \item{head}{\code{signature(object = "phylo4")}: show first few
      nodes}
    \item{tail}{\code{signature(object = "phylo4")}: show last few
      nodes}  
  }
}

\usage{
    \S4method{print}{phylo4}(x, edgeOrder = c("pretty", "real"), printall)
    \S4method{summary}{phylo4}(object, quiet=FALSE)
    \S4method{show}{phylo4}(object)
    \S4method{names}{phylo4}(x)
}

\arguments{
    \item{x}{a phylo4 object}
    \item{object}{a phylo4 object}
    \item{edgeOrder}{Character string indicating whether the edges
      should be printed as ordered in the tree "real" (e.g. preorder or
      postorder), or "pretty" printed with tips collated together}
    \item{printall}{If TRUE all tip labels are printed}
    \item{quiet}{a logical stating whether the results of the summary
      should be printed to the screen (FALSE, default) or not (TRUE)}
}
  
\value{
    
  The \code{summary} method invisibly returns a list with the
  following components:
  
  \item{\code{name}}{the name of the object}
  \item{\code{nb.tips}}{the number of tips}
  \item{\code{nb.nodes}}{the number of nodes}
  \item{\code{mean.el}}{mean of edge lengths}
  \item{\code{var.el}}{variance of edge lengths (estimate for population)}
  \item{\code{sumry.el}}{summary (i.e. range and quartiles) of the edge
    lengths}
  \item{\code{degree}}{(optional) degree (i.e. number of descendants) of
    each node; displayed only when there are polytomies}
  \item{\code{polytomy}}{(optional) type of polytomy for each node:
    \sQuote{node}, \sQuote{terminal} (all descendants are tips) or
    \sQuote{internal} (at least one descendant is an internal node);
    displayed only when there are polytomies}
  
  The \code{names} method returns a vector of characters corresponding
  to the names of the slots.
}

\seealso{
  The \code{\link{phylo4}} constructor, the \code{\link{checkPhylo4}} 
  function to check the validity of \code{phylo4} objects. See also the
  \code{\link{phylo4d}} constructor and the \linkS4class{phylo4d}
  class.
} 
\author{Ben Bolker, Thibaut Jombart}
\examples{

  tOwls <- "(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
  tree.owls <- read.tree(text=tOwls)
  P1 <- as(tree.owls, "phylo4")
  P1
  summary(P1)


  ## summary of a polytomous tree
  E <- matrix(c(
      8,  9,
      9, 10,
     10,  1,
     10,  2,
      9,  3,
      9,  4,
      8, 11,
     11,  5,
     11,  6,
     11,  7,
      0,  8), ncol=2, byrow=TRUE)

  P2 <- phylo4(E)
  nodeLabels(P2) <- as.character(nodeId(P2, "internal"))
  plot(P2, show.node.label=TRUE)
  sumryP2 <- summary(P2)
  sumryP2
}
\keyword{methods}
