\name{readNexus}
\docType{methods}
\alias{readNexus}
\title{Parse a Nexus file for one or more tree(s) and/or tip data}
\description{
  \code{readNexus} reads a Nexus file and outputs a \code{data.frame}, \code{phylo4}, or \code{phylo4d} object.
}
\section{Methods}{
  \describe{
\item{x = "readNexus"}{creates a \code{data.frame}, \code{phylo4}, or \code{phylo4d} object from a Nexus file}

   }
 }
 \usage{
readNexus(file, simplify=TRUE, type=c("all","tree","data"), char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=TRUE)
}
\arguments{
  \item{file}{a Nexus file}
  \item{simplify}{If FALSE, returns a list. If TRUE, returns a list if there are multiple trees, a single object otherwise}
  \item{type}{Determines which type of objects to return, if present}
  \item{char.all}{If TRUE, returns all characters, even those excluded in the NEXUS file}
  \item{polymorphic.convert}{If TRUE, converts polymorphic DNA characters to missing data}
  \item{levels.uniform}{If TRUE, uses the same levels for all characters}
}
\value{
  A \linkS4class{phylo4} tree, \linkS4class{phylo4d} object, or a list of trees if simplify=FALSE or there is more than one tree, or a data.frame.
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
