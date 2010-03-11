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
	readNexus(file, simplify=FALSE, type=c("all", "tree", "data"),
char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=TRUE,
quiet=TRUE, check.node.labels=c("keep", "drop", "asdata"), return.labels=TRUE, ...)
 }
\arguments{
  \item{file}{a Nexus file}
  \item{simplify}{If FALSE, returns a list. If TRUE, returns a list if
    there are multiple trees, a single object otherwise} 
  \item{type}{Determines which type of objects to return, if present}
  \item{char.all}{If TRUE, returns all characters, even those excluded
    in the NEXUS file} 
  \item{polymorphic.convert}{If TRUE, converts polymorphic DNA
    characters to missing data} 
  \item{levels.uniform}{If TRUE, uses the same levels for all
    characters}
  \item{quiet}{If FALSE the tree string is printed. This is mainly for
    debugging purposes. This option can considerably slow down the
    process if the tree is big or there are many trees in the file.}
  \item{check.node.labels}{Determines how the node labels in the Nexus
    files should be treated in the phylo4 object, see Details for more
    information.}
  \item{return.labels}{Determines whether state names (if TRUE) or state
    codes should be returned.}
  \item{\dots}{Additional arguments to be passed to phylo4 or phylo4d
    constructor (see details)}
}
\details{
  The options for \code{check.node.labels} can take the values:
  \describe{
     \item{keep}{the node labels of the Nexus file will be passed as
       node labels in the phylo4 object} 
     \item{drop}{the node labels of the Nexus file will be ignored in
       the phylo4 object} 
     \item{asdata}{the node labels will passed to the phylo4 object as
       data.} 
 }
 If you use the option \code{asdata} on a file with no node labels, a
 warning message is issued and by default \code{check.node.labels} takes
 the value \code{drop}. 

 Depending on the context \code{readNexus} will call either the
 \code{phylo4} or \code{phylo4d} constructor. The \code{phylo4d}
 constructor will be used if \code{check.node.labels="asdata"} or if
 \code{type="all"}; and the \code{phylo4} constructor will be used
 otherwise. Additional arguments can be passed to the constructors such
 as \code{annote}, \code{missing.data} or \code{extra.data}. See the
 documentation of \link{phylo4-methods}, \link{phylo4d} and
 \link{formatData} for the complete list of arguments.
 
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
  This relies on Version 2.0 of the Nexus Class Library by Paul Lewis
  and Mark Holder. Tree reading is done initially using internal
  \code{APE} functions; these functions can be confused if tree names
  include the word (\dQuote{tree}) in them.
}
\keyword{misc}
