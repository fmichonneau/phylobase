\docType{class}
\name{phylo4d-class}
\alias{phylo4d-class}
\title{phylo4d class}
\description{
S4 class for phylogenetic tree and data.
}
\section{Objects from the Class}{
  Objects can be created from various trees and a
  data.frame using the constructor \code{phylo4d}, or using
  \code{new("phylo4d", \dots{})} for empty objects.
}
\examples{
example(read.tree, "ape")
  obj <- phylo4d(as(tree.owls.bis,"phylo4"), data.frame(wing=1:3))
  obj
  names(obj)
  summary(obj)
}
\author{
Ben Bolker, Thibaut Jombart
}
\seealso{
\code{\link{coerce-methods}} for translation functions. The
\code{\link{phylo4d}} constructor and the
\code{\link{formatData}} function to check the validity of
trees and data. See also the \code{\link{phylo4}}
constructor, the \linkS4class{phylo4} class, and the
\code{\link{checkPhylo4}} function to check the validity of
\code{phylo4} trees.
}
\keyword{classes}

