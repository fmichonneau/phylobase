\name{checkPhylo4}
\alias{checkPhylo4}
\alias{checkTree}
\alias{checkPhylo4Data}
\title{Validity checking for phylo4 objects}
\description{
  Basic checks on the validity of S4 phylogenetic objects
}
\usage{
checkPhylo4(object)
checkTree(object,
             warn=c("retic","singleton","multiroot"),
             err=NULL)
checkPhylo4Data(object)
}
\arguments{
  \item{object}{A prospective phylo4 or phylo4d object}
  \item{warn}{a character vector listing phenomena to warn about:
    current options are [\code{"poly"} (polytomies),
    \code{"retic"} (reticulations),
    \code{"singleton"} ("singleton" nodes, i.e.
    nodes with a single descendant), \code{"multiroot"}
    (multiple roots) ]}
  \item{err}{a character vector listing phenomena to trigger errors:
    choices as above}
}
\value{
  As required by \code{\link[methods]{validObject}}, returns an
  error string (describing problems) or TRUE if everything is OK.
}
\note{

  These functions are only intended to be called by other phylobase
  functions.

  \code{checkPhylo4} is an (inflexible) wrapper for
  \code{checkTree}, which uses the default settings.
  The rules for \code{phylo4} objects essentially follow
  those for \code{phylo} objects from the \code{ape} package,
  which are in turn defined in
  http://ape.mpl.ird.fr/misc/FormatTreeR_4Dec2006.pdf.
  These are essentially that:
  \itemize{
    \item if the tree has edge lengths defined, the number of edge
    lengths must match the number of edges; 
    \item the number of tip labels must match the number of tips;
    \item in a tree with \code{ntips} tips and \code{nnodes} (total)
    nodes, nodes 1 to \code{ntips} must be tips 
    \item if the tree is rooted, the root must be node number
    \code{ntips+1} and the root node must be the first row of the edge
    matrix
    \item tip labels, node labels, edge labels, edge lengths must have
    proper internal names (i.e. internal names that match the node
    numbers they document)
    \item tip and node labels must be unique
  }

  For \code{phylo4d} objects, \code{checkTree} also calls
  \code{checkPhylo4Data} to check the validity of the data associated
  with the tree. It ensures that (1) the data associated with the tree
  have the correct dimensions, (2) that the row names for the data are
  correct.
}

\seealso{
  the \code{\link{phylo4}} constructor and \linkS4class{phylo4} class;
  \code{\link{formatData}},  the \code{\link{phylo4d}} constructor and
  the \linkS4class{phylo4d} class do checks for the data associated with trees.
  See \code{\link{coerce-methods}} for translation functions.
}
\author{Ben Bolker, Steven Kembel, Francois Michonneau}
\keyword{misc}
