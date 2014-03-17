\name{phylomat-class}
\docType{class}
\alias{phylo4vcov-class}
\alias{as_phylo4vcov}
\title{matrix classes for phylobase}
\description{Classes representing phylogenies as matrices}
\section{Objects from the Class}{
  These are square matrices (with rows and columns corresponding
  to tips, and internal nodes implicit) with different
  meanings depending on the type (variance-covariance matrix,
  distance matrix, etc.).
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{square, numeric matrix with row and column
      labels corresponding to the tip labels}
    \item{\code{edge.label}:}{character vector of edge labels}
    \item{\code{order}:}{character describing original ordering of edge matrix}
  }
}
\usage{
as_phylo4vcov(from, \dots)
}
\arguments{
  \item{from}{a \code{phylo4} object}
  \item{\dots}{optional arguments, to be passed to \code{vcov.phylo} in
    \code{ape} (the main useful option is \code{cor}, which can be set
    to \code{TRUE} to compute a correlation rather than a
    variance-covariance matrix)}
}
\author{Ben Bolker}
\examples{
   trString <- "(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);"
   tree.owls <- ape::read.tree(text=trString)
  o2 <- as(tree.owls,"phylo4")
  ov <- as(o2,"phylo4vcov")
  o3 <- as(ov,"phylo4")
  ## these are not completely identical, but are
  ## topologically identical ...

  ## edge matrices are in a different order:
  ## cf. o2@edge and o3@edge
  ## BUT the edge matrices are otherwise identical
  identical(o2@edge[order(o2@edge[,2]),],
            o3@edge[order(o3@edge[,2]),])

  ## There is left/right ambiguity here in the tree orders:
  ## in o2 the 5->6->7->1 lineage
  ## (terminating in Strix aluco)
  ## is first, in o3 the 5->6->3 lineage
  ## (terminating in Athene noctua) is first.

}
\keyword{classes}
