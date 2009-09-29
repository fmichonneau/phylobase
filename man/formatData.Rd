\name{formatData}
\alias{formatData}

\title{Formatting data before attaching them to a phylo4 object}

\description{
  Basic checks on the validity of S4 phylogenetic plus data objects
  (i.e. \linkS4class{phylo4d} objects).
}

\usage{
  formatData(phy, dt, type=c("tip", "internal", "all"),
              match.data=TRUE, rownamesAsLabels=FALSE,
              label.type=c("rownames", "column"),
              label.column=1, missing.data=c("fail", "warn", "OK"),
              extra.data=c("warn", "OK", "fail"), keep.all=TRUE)
}

\arguments{
  \item{phy}{a valid \code{phylo4} object}
  \item{dt}{a data frame or a vector}
  \item{type}{type of data to attach}
  \item{match.data}{(logical) Are data labels matched against tree
    labels?} 
  \item{label.type}{character, \code{rownames} or \code{column}: should
    the labels be taken from the row names of \code{dt} or from
    the \code{label.column} column of \code{dt}?}
  \item{label.column}{if \code{label.type=="column"}, column specifier
    (number or name) of the column containing tip labels}
  \item{missing.data}{action to take if there are missing data or if
    there are data labels that don't match}
  \item{extra.data}{action to take if there are extra data or if there
    are labels that don't match} 
  \item{rownamesAsLabels}{(logical), in the case of number-like labels
    should the row names of the data provided be considered as labels
    (TRUE) or node numbers (FALSE and default)}
  \item{keep.all}{(logical), should the data have rows for all nodes
    (with NA values for internal rows when type='tip', and vice versa) 
    (TRUE and default) or only rows corresponding to type argument}
}

\value{
  \code{formatData} returns a data frame having node numbers as row
  names. The data frame is also formatted to have the correct dimension
  given the \code{phylo4} object provided.
}

\details{

  \code{formatData} is an internal function that should not be called by
  the user. It formats properly the data provided by the user before
  being attached to a tree.
  
  By default, data names are used to attach tip data to the
  tree. In this case data names must be present and must match the set
  of tip or internal nodes exactly. The matching process can be done
  either by using \code{phylo4} labels (i.e. tip or node) or by using
  node numbers. It is possible to provide both at the same time by using
  the option \code{rownamesAsLabels=FALSE}.
  
  \code{formatData} (1) converts labels provided in the data into node
  numbers, (2) makes sure that the data are appropriately matched
  against tip or nodes, (3) checks for differences between data and
  tree, (4) creates a data frame with the correct dimensions given a
  tree.
  
}

\seealso{the \code{\link{phylo4d}} constructor, the
  \linkS4class{phylo4d} class. See also the \code{\link{checkPhylo4}},
  the \code{\link{phylo4}} constructor and the \linkS4class{phylo4}
  class. See \code{\link{coerce-methods}} for translation functions.}

\author{Francois Michonneau}

\keyword{misc}
