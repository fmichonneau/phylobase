\name{readNexus}
\docType{methods}
\alias{readNexus}
\title{Create a phylo4 or a phylo4d object from a Nexus file}
\description{
  \code{readNexus} reads a Nexus file and outputs a \code{phylo4} or \code{phylo4d} object.
}
\section{Methods}{
  \describe{
	\item{x = "readNexus"}{creates a \linkS4class{phylo4} or \linkS4class{phylo4d} object from a Nexus file}
   }
 }
 \usage{
	readNexus(file, simplify=TRUE, type=c("all","tree","data"),
		   char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=TRUE,
                   check.node.labels=c("keep", "drop", "asdata"))
 }
\arguments{
  \item{file}{a Nexus file}
  \item{simplify}{If FALSE, returns a list. If TRUE, returns a list if there are multiple trees, a single object otherwise}
  \item{type}{Determines which type of objects to return, if present}
  \item{char.all}{If TRUE, returns all characters, even those excluded in the NEXUS file}
  \item{polymorphic.convert}{If TRUE, converts polymorphic DNA characters to missing data}
  \item{levels.uniform}{If TRUE, uses the same levels for all characters}
  \item{check.node.labels}{Determines how the node labels in the Nexus files should be
 treated in the phylo4 object, see Details for more information.}
}
\details{
  The options for \code{check.node.labels} can take the values:
  \describe{
     \item{keep}{the node labels of the Nexus file will be passed as node labels in the phylo4 object}
     \item{drop}{the node labels of the Nexus file will be ignored in the phylo4 object}
     \item{asdata}{the node labels will passed to the phylo4 object as data.}
 }
 If you use the option \code{asdata} on a file with no node labels, a warning message is
 issued and by default \code{check.node.labels} takes the value \code{drop}.
}
\value{
  A \linkS4class{phylo4} object if there is no data in the Nexus file
  and if the option \code{check.node.labels} is \dQuote{keep} or
  \dQuote{drop}, or a \linkS4class{phylo4d} object otherwise.
  If several trees are included in the Nexus file and the option
  \code{simplify=FALSE} a list of \linkS4class{phylo4} or
  \linkS4class{phylo4d} objects is returned.
}
\seealso{the \linkS4class{phylo4d} class, the \linkS4class{phylo4} class}
\author{Brian O'Meara, Derrick Zwickl}
\note{
  This relies on Version 2.0 of the Nexus Class Library by Paul Lewis and Mark Holder. Tree reading
  is done initially using internal \code{APE} functions; these functions can be
  confused if tree names include the word (\dQuote{tree}) in them. 
}
\keyword{misc}
