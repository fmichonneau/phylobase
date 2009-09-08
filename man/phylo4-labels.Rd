\name{labels-phylo4}
\docType{methods}
\alias{labels<-}
\alias{labels,phylo4-method}
\alias{labels<-,phylo4,ANY,ANY,character-method}
\alias{labels<-,phylo4d,ANY,ANY,character-method}
\alias{hasNodeLabels}
\alias{hasNodeLabels-methods}
\alias{hasNodeLabels,phylo4-method}
\alias{nodeLabels}
\alias{nodeLabels-methods}
\alias{nodeLabels,phylo4-method}
\alias{nodeLabels<-}
\alias{nodeLabels<-,phylo4,character-method}
\alias{nodeLabels<-,phylo4d,ANY-method}
\alias{tipLabels}
\alias{tipLabels-methods}
\alias{tipLabels,phylo4-method}
\alias{tipLabels<-}
\alias{tipLabels<-,phylo4,character-method}
\alias{tipLabels<-,phylo4d,character-method}
\alias{hasEdgeLabels}
\alias{hasEdgeLabels-methods}
\alias{hasEdgeLabels,phylo4-method}
\alias{edgeLabels}
\alias{edgeLabels<-}
\alias{edgeLabels-methods}
\alias{edgeLabels,phylo4-method}
\alias{edgeLabels<-,phylo4,character-method}

\title{Labels for phylo4/phylo4d objects}

\description{Methods for creating, accessing and updating labels in
  phylo4/phylo4d objects}

\usage{

\S4method{labels}{phylo4}(object, type=c("all", "tip", "internal"), \dots)
\S4method{labels}{phylo4,ANY,ANY,character}(x, type=c("tip", "internal", "allnode"),
use.names=FALSE) <- value

\S4method{tipLabels}{phylo4}(object)
\S4method{tipLabels}{phylo4,character}(object) <- value

\S4method{hasNodeLabels}{phylo4}(x)
\S4method{nodeLabels}{phylo4}(object)
\S4method{nodeLabels}{phylo4,character}(object) <- value

\S4method{hasEdgeLabels}{phylo4}(x)
\S4method{edgeLabels}{phylo4}(x)
\S4method{edgeLabels}{phylo4,character}(x) <- value

}

\arguments{
  \item{x}{a phylo4 or phylo4d object.}
  \item{object}{a phylo4 or phylo4d object.}
  \item{type}{which type of labels: \code{all} (tips and internal nodes),
    \code{tip} (tips only), \code{internal} (internal nodes only).}
  \item{value}{a vector of class \code{character}, see Details for more
    information.}
  \item{use.names}{should the names of the vector used to create/update
    labels be used to match the labels? See Details for more information.}
  \item{\dots}{additional parameter passed to the function (currently ignored)}
}

\details{

  In phylo4/phylo4d objects, tips must have labels (that's why there is
  no method for hasTipLabels), internal nodes and edges can have
  labels. 

  Labels must be provided as a vector of class \code{character}. The
  length of the vector must match the number of elements they label.

  The option \code{use.names} allows the user to match a label to a
  particular node. In this case, the vector must have names that match
  the node numbers.

  The function \code{labels} is mostly intended to be used
  internally.  

}

\section{Methods}{
  \describe{
    \item{labels}{\code{signature(object = "phylo4")}: tip and/or
      internal node labels, ordered by node ID}

    \item{tipLabels}{\code{signature(object = "phylo4")}: tip labels,
      ordered by node ID}
    
    \item{hasNodeLabels}{\code{signature(object = "phylo4")}: whether
      tree has (internal) node labels}
    \item{nodeLabels}{\code{signature(object = "phylo4")}: internal
      node labels, ordered by node ID}
    
    \item{hasEdgeLabels}{\code{signature(object = "phylo4")}: whether
      tree has (internal) edge labels}
    \item{edgeLabels}{\code{signature(object = "phylo4")}: internal
      edge labels, ordered according to the edge matrix}
  }
}

\examples{

data(geospiza)

## Return labels from geospiza
tipLabels(geospiza)

## Internal node labels in geospiza are empty
nodeLabels(geospiza)

## Creating internal node labels
ndLbl <- paste("n", 1:nNodes(geospiza), sep="")
nodeLabels(geospiza) <- ndLbl
nodeLabels(geospiza)

## naming the labels
names(ndLbl) <- nodeId(geospiza, "internal")

## shuffling the labels
(ndLbl <- sample(ndLbl))

## by default, the labels are attributed in the order
## they are given:
nodeLabels(geospiza) <- ndLbl
nodeLabels(geospiza)

## but use.names puts them in the correct order
labels(geospiza, "internal", use.names=TRUE) <- ndLbl
nodeLabels(geospiza)

}
