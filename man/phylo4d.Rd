\name{phylo4d-methods}
\docType{methods}
\alias{phylo4d}
\alias{phylo4d-methods}
\alias{phylo4d,phylo4-method}
\alias{phylo4d,matrix-method}
\alias{phylo4d,phylo-method}
\title{Combine a phylogenetic tree with data}
\description{
  \code{phylo4d} is a generic constructor which merges a phylogenetic tree with
  data frames to create a combined object of class \code{phylo4d}
}
\section{Methods}{
  \describe{
     \item{x = "phylo4"}{merges a tree of class \code{phylo4} with a
       data.frame into a \code{phylo4d} object}

     \item{x = "matrix"}{merges a matrix of tree edges similar to the edge
       slot of a \code{phylo4} object (or to \$edge of a \code{phylo} object) with a
       data.frame into a \code{phylo4d} object}

     \item{x = "phylo"}{merges a tree of class \code{phylo} with a
       data.frame into a \code{phylo4d} object }    
   }}
 
\usage{
     \S4method{phylo4d}{phylo4}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL, merge.tip.node = TRUE, ...)
     \S4method{phylo4d}{phylo}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL, merge.tip.node = TRUE, ...)
     \S4method{phylo4d}{matrix}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL, merge.tip.node = TRUE, ...)
}

\arguments{
  \item{x}{an object of class \code{phylo4}, \code{phylo} or a matrix of
  edges (see above)}
  \item{tip.data}{a data frame for tips data}
  \item{node.data}{a data frame for nodes data}
  \item{all.data}{a data frame for all (i.e. tips and nodes) data. In such
    case, first rows should correspond to tips, last rows to nodes.}
  \item{merge.tip.node}{if both \code{tip.data} and \code{node.data} are
    provided, it determines if they should be merged as a single
    trait. This argument is evaluated only if both \code{tip.data} and
    \code{node.data} have identical column names.}
  \item{\dots}{further arguments to be passed to other methods. Used to
    control the validity checking of the created object (see \code{\link{check_data}}).}
}

\details{
  You can provide several data frames to define traits associated with
  tips and/or nodes. If you provide \code{all.data} and \code{tip.data}
  or \code{node.data}, row names of the data frames will be matched
  (\code{all.data} names are matched against \code{tip.data} and/or
  \code{node.data}). This is done independently of the labels of the tree 
  (and also of the value of the arguments \code{use.tip.names} and
  \code{use.node.names}). This means that you need to be consistent
  with the row names of your data frames. It is good practice to use tip
  and node labels when you merge data with a tree.
  If \code{all.data} is provided and tip and node names are not being
  used, then the first part of the data frame will be tips data and the
  second part will be nodes data.
}

\value{
  An object of class \linkS4class{phylo4d}.
}

\seealso{
\code{\link{coerce-methods}} for translation functions. The
  \linkS4class{phylo4d} class, the \code{\link{check_data}} 
  function to check the validity of \code{phylo4d} objects;
  \linkS4class{phylo4} class and \link{phylo4} constructor.}

\author{Ben Bolker, Thibaut Jombart, Steve Kembel, \enc{François}{Francois} Michonneau}

\note{
  Checking on matches between the tree and the data will be done by the
  validity checker (label matches between data and tree tips, number of
  rows of data vs. number of nodes/tips/etc.)
}
\examples{
library(ape)
example(read.tree)
try(phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3)), silent=TRUE)
obj <- phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3),
use.tip.names=FALSE)
obj
print(obj)

####

data(geospiza_raw)
geoTree <- geospiza_raw$tree
geoData <- geospiza_raw$data

## fix differences in tip names between the tree and the data
geoData <- rbind(geoData, array(, dim = c(1,ncol(geoData)),
                  dimnames = list("olivacea", colnames(geoData))))

### Example using a tree of class 'phylo'
exGeo1 <- phylo4d(geoTree, tip.data = geoData)

### Example using a tree of class 'phylo4'
geoTree <- as(geoTree, "phylo4")

## some random node data
rNodeData <- data.frame(randomTrait = rnorm(nNodes(geoTree)))
## comment out node labels for now (empty)
##                        row.names = nodeLabels(geoTree))

exGeo2 <- phylo4d(geoTree, tip.data = geoData, node.data = rNodeData)

### Example using 'merge.tip.node'
## some random tip data
rTipData <- data.frame(randomTrait = rnorm(nTips(geoTree)),
                       row.names = labels(geoTree))

(exGeo3 <- phylo4d(geoTree, tip.data = rTipData, node.data = rNodeData))
(exGeo4 <- phylo4d(geoTree, tip.data = rTipData, node.data = rNodeData,
                    merge.tip.node = FALSE))

		    ### Example with 'all.data'
nodeLabels(geoTree) <- as.character(nodeId(geoTree))
rAllData <- data.frame(randomTrait = rnorm(nTips(geoTree) + nNodes(geoTree)),
row.names = c(labels(geoTree),nodeId(geoTree)))


exGeo5 <- phylo4d(geoTree, all.data = rAllData)
       

}
\keyword{misc}
