\name{NexusToDataFrame}
\docType{methods}
\alias{NexusToDataFrame}
\title{Load data from a Nexus file}
\description{
  \code{NexusToDataFrame} reads a Nexus file and outputs a data.frame with data.
}
\section{Methods}{
  \describe{
\item{x = "NexusToDataFrame"}{creates a data.frame from a Nexus file}

   }
 }
 \usage{
NexusToDataFrame(fileToRead, allchar = FALSE, polymorphictomissing = TRUE, levelsall = TRUE)
}
\arguments{
  \item{fileToRead}{a Nexus file}
\item{allchar}{if TRUE, load all characters, even excluded ones}
\item{polymorphictomissing}{if TRUE, convert polymorphic data to missing data}
\item{levelsall}{categorical data are loaded with the same levels, even
  if one character is missing a state}
}
\value{
  A dataframe suitable to use with \linkS4class{phylo4d}.
}
\seealso{the \linkS4class{phylo4d} class, the \linkS4class{phylo4} class}
\author{Brian O'Meara, Derrick Zwickl}
\note{
  This relies on Version 2.0 of the Nexus Class Library by Paul Lewis and Mark Holder.
  If using Nexus files from Mesquite with continuous data, there must not be an 
  (\dQuote{Items = ()}) command in the file; delete it if present.
}
\keyword{misc}
