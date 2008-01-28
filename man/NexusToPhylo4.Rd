\name{NexusToPhylo4}
\docType{methods}
\alias{NexusToPhylo4}
\title{Load a tree from a Nexus file}
\description{
  \code{NexusToPhylo4} reads a Nexus file and outputs a \code{phylo4} tree.
}
\section{Methods}{
  \describe{
\item{x = "NexusToPhylo4"}{creates a \code{phylo4} tree from a Nexus file}

   }
 }
 \usage{
\NexusToPhylo4(fileToRead,multi=FALSE)
}
\arguments{
  \item{fileToRead}{a Nexus file}
\item{multi}{if TRUE, return a list of \code{phylo4} trees}
}
\value{
  A \linkS4class{phylo4} tree, or a list of trees if multi=TRUE.
}
\seealso{the \linkS4class{phylo4d} class, the \linkS4class{phylo4} class}
\author{Brian O'Meara, Derrick Zwickl}
\note{
  This relies on Version 2.0 of the Nexus Class Library by Paul Lewis and Mark Holder.
  Tree reading is done initially using a function from APE; this function can be 
  confused if the tree name has the word (\dQuote{tree}) in it. Also note that edge
  and node labels (except for terminal nodes) are NOT yet loaded into the \code{phylo4}
  object.
}
\keyword{misc}
