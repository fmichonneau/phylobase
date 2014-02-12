\name{Import Nexus and Newick files}
\docType{methods}
\alias{readNexus}
\alias{readNCL}
\alias{readNewick}

\title{Create a phylo4, phylo4d or data.frame object from a Nexus or a
  Newick file}

\description{
  \code{readNexus} reads a Nexus file and outputs a \code{phylo4} or
  \code{phylo4d} or \code{data.frame} object. 
  
  \code{readNewick} reads a Newick file and outputs a \code{phylo4} or
  \code{phylo4d} object. 
}

\usage{
readNexus(file, simplify=FALSE, type=c("all", "tree", "data"),
          char.all=FALSE, polymorphic.convert=TRUE, levels.uniform=FALSE,
          quiet=TRUE, check.node.labels=c("keep", "drop", "asdata"),
          return.labels=TRUE, check.names=TRUE, convert.edge.length=FALSE,
          ...)

readNewick(file, simplify=FALSE, quiet=TRUE,
           check.node.labels=c("keep", "drop", "asdata"),
           convert.edge.length=FALSE, ...)

}

\arguments{
  \item{file}{a Nexus file for \code{readNexus} or a file that contains
  Newick formatted trees for \code{readNewick}}
  
  \item{simplify}{If there are multiple trees in the file, only the
    first one is returned if TRUE and a list of phylo4/phylo4d objects
    is returned if the file contains multiple trees.}
  
  \item{type}{Determines which type of objects to return, if present in 
    the file (see Details).}
  
  \item{char.all}{If TRUE, returns all characters, even those excluded
    in the NEXUS file}
  
  \item{polymorphic.convert}{If TRUE, converts polymorphic
    characters to missing data}
  
  \item{levels.uniform}{If TRUE, uses the same levels for all
    characters}
  
  \item{quiet}{If FALSE the output of the NCL interface is printed. This
    is mainly for debugging purposes. This option can considerably slow
    down the  process if the tree is big or there are many trees in the
    file.}
  
  \item{check.node.labels}{Determines how the node labels in the Nexus
    or Newick files should be treated in the phylo4 object, see Details
    for more information.}
  
  \item{return.labels}{Determines whether state names (if TRUE) or state
    codes should be returned.}

  \item{check.names}{logical. If \sQuote{TRUE} then the names of the
          characters from the NEXUS file are checked to ensure that they
          are syntactically valid variable names and are not duplicated.
          If necessary they are adjusted (by \sQuote{make.names}) so
          that they are.}

  \item{convert.edge.length}{logical. If \sQuote{TRUE} negative edge
    lengths are replaced with 0. At this time \code{phylobase} does not
    accept objects with negative branch lengths, this workaround allows
    to still use trees with negative branch lengths are an artifact of
    the method used to build the tree.}
  
  \item{\dots}{Additional arguments to be passed to phylo4 or phylo4d
    constructor (see Details)}
}
\details{

 \code{readNexus} extracts data held in a Nexus file, specifically from
 DATA, CHARACTER or TREES blocks present in the file. The \code{type} 
 argument specifies which of these is returned:
 \describe{
 \item{data}{will only return a \code{data.frame} of the contents of all 
   DATA and CHARACTER blocks.}
 \item{tree}{will only return a \code{phylo4} object of the contents of 
   the TREES block.}
 \item{all}{if only data or a tree are present in the file, this option 
   will act as the options above, returning either a \code{data.frame} 
   or a \code{phylo4} object respectively. If both are present then a
   \code{phylo4d} object is returned containing both.}
 }
 The function returns \code{NULL} if the \code{type} of data requested is 
 not present in the file, or if neither data nor tree blocks are present.

 Depending on the context \code{readNexus} will call either the
 \code{phylo4} or \code{phylo4d} constructor. In addition with
 \code{type="all"}, the \code{phylo4d} constructor will be used if 
 \code{check.node.labels="asdata"}. 

 \code{readNewick} imports newick formatted tree files and will return a
 \code{phylo4} or a \code{phylo4d} object if the option
 \code{check.node.labels="asdata"} is invoked.
 
 For both \code{readNexus} and \code{readNewick}, the options for
 \code{check.node.labels} can take the values:
 \describe{
   \item{keep}{the node labels of the trees will be passed as node
     labels in the \code{phylo4} object}
   \item{drop}{the node labels of the trees will be ignored in
     the \code{phylo4} object} 
   \item{asdata}{the node labels will be passed as data and a
     \code{phylo4d} object will be returned.} 
 }
 If you use the option \code{asdata} on a file with no node labels, a
 warning message is issued, and thus \code{check.node.labels} takes
 the value \code{drop}. 

 For both \code{readNexus} and \code{readNewick}, additional arguments
 can be passed to the constructors such as \code{annote},
 \code{missing.data} or \code{extra.data}. See the documentation of
 \link{phylo4-methods}, \link{phylo4d} and \link{formatData} for the
 complete list of options.
 
}
\note{Underscores in state labels (i.e. trait or taxon names)
  will be translated to spaces when read by NCL. Unless
  \code{check.names=FALSE}, trait names will
  be converted to valid R names (see \code{\link{make.names}}) on
  input to R, so spaces will be translated to periods.}
\value{
  Depending on the value of \code{type} and the contents of the file, 
  one of: a \code{data.frame}, a \linkS4class{phylo4} object,  a
  \linkS4class{phylo4d} object or \code{NULL}.
  If several trees are included in the Nexus file and the option
  \code{simplify=FALSE} a list of \linkS4class{phylo4} or
  \linkS4class{phylo4d} objects is returned.
}
\seealso{the \linkS4class{phylo4d} class, the \linkS4class{phylo4} class}
\author{Brian O'Meara, Francois Michonneau, Derrick Zwickl}
\note{
  This relies on Version 2.1.13 of the Nexus Class Library by Paul Lewis
  and Mark Holder. Tree reading is done initially using internal
  \code{APE} functions; these functions can be confused if tree names
  include the word (\dQuote{tree}) in them.
}
\keyword{misc}
