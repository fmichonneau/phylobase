\name{NexusToPhylo4D}
\docType{methods}
\alias{NexusToPhylo4D}
\title{Load data and a tree from a Nexus file}
\description{
  \code{NexusToPhylo4D} reads a Nexus file and outputs a \code{phylo4d} object.
}
\section{Methods}{
  \describe{
\item{x = "NexusToPhylo4D"}{creates a \code{phylo4d} object from a Nexus file}

   }
 }
 \usage{
\NexusToPhylo4D(fileToRead, multi=FALSE, allchar=FALSE, polymorphictomissing=TRUE, levelsall=TRUE)
}
\arguments{
  \item{fileToRead}{a Nexus file}
\item{multi}{if TRUE, return a list of \code{phylo4} trees}
\item{allchar}{if TRUE, load all characters, even excluded ones}
\item{polymorphictomissing}{if TRUE, convert polymorphic data to missing data}
\item{levelsall}{categorical data are loaded with the same levels, even
  if one character is missing a state}
}
\value{
  A \linkS4class{phylo4d} object.
}
\seealso{the \linkS4class{phylo4d} class, the \linkS4class{phylo4} class}
\author{Brian O'Meara, Derrick Zwickl}
\note{
  This relies on Version 2.0 of the Nexus Class Library by Paul Lewis and Mark Holder.
  If using Nexus files from Mesquite with continuous data, there must not be an 
  (\dQuote{Items = ()}) command in the file; delete it if present.
  Tree reading is done initially using a function from APE; this function can be 
  confused if the tree name has the word (\dQuote{tree}) in it. Also note that edge
  and node labels (except for terminal nodes) are NOT yet loaded into the \code{phylo4d}
  object.
}
\keyword{misc}
