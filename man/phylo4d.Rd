\name{phylo4d}
\docType{methods}
\alias{phylo4d}
\alias{phylo4d-methods}
\alias{phylo4d,phylo4-method}
\alias{phylo4d,phylo4d-method}
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
     \S4method{phylo4d}{phylo}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL, check.node.labels = c("keep", "drop", "asdata"), 
        annote=list(), metadata=list(), ...)
     \S4method{phylo4d}{phylo4}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL, match.data = TRUE, merge.data = TRUE, 
        rownamesAsLabels = FALSE, metadata = list(), ...)
     \S4method{phylo4d}{matrix}(x, tip.data = NULL, node.data = NULL,
        all.data = NULL, metadata = list(), ...)
}

\arguments{
  \item{x}{an object of class \code{phylo4}, \code{phylo} or a matrix of
  edges (see above)}
  \item{tip.data}{a data frame (or object to be coerced to one)
    containing only tip data}
  \item{node.data}{a data frame (or object to be coerced to one)
    containing only node data}
  \item{all.data}{a data frame (or object to be coerced to one)
    containing both tip and node data}
  \item{match.data}{(logical) should the rownames of the data frame
    provided in \code{tip.data}, \code{node.data} and/or \code{all.data}
    be used to be matched against tip and internal node identifiers? See
    details.}
  \item{merge.data}{if both \code{tip.data} and \code{node.data} are
    provided, should columns with common names will be merged together
    (default TRUE) or not (FALSE)? See details.}
  \item{check.node.labels}{if \code{x} is of class \code{phylo}, use
    either "keep" (the default) to retain internal node labels, "drop"
    to drop them, or "asdata" to convert them to numeric tree data. This
    argument is useful if the \code{phylo} object has non-unique node
    labels or node labels with informative data (e.g., posterior
    probabilities).}
  \item{rownamesAsLabels}{(logical), should the row names of the data
    provided be matched only to labels (TRUE), or should any number-like
    row names be matched to node numbers (FALSE and default)}
  \item{annote}{any additional annotation data to be passed to the new object}
  \item{metadata}{any additional metadata to be passed to the new object}
  \item{\dots}{further arguments to be passed to
    \code{\link{formatData}}. Notably, these additional arguments
    control the behavior of the constructor in the case of missing/extra
    data and where to look for labels in the case of non-unique labels
    that cannot be stored as row names in a data frame.}
}

\details{
  You can provide several data frames to define traits associated with
  tip and/or internal nodes. Rules for matching rows of data to tree
  nodes are determined jointly by the \code{match.data} and
  \code{rownamesAsLabels} arguments, and these will apply the same way
  to all supplied data frames. If \code{match.data} is TRUE, data frame
  rows will be matched exclusively against tip and node labels if
  \code{rownamesAsLabels} is also TRUE, whereas any all-digit row names
  will be matched against tip and node numbers if
  \code{rownamesAsLabels} is FALSE (the default). If \code{match.data}
  is FALSE, \code{rownamesAsLabels} has no effect, and row matching is
  purely positional with respect to the order returned by
  \code{nodeId(x, "tip")} for tip.data, \code{nodeId(x, "internal")} for
  node.data, and \code{nodeId(x, "all")} for all.data.

  This means that you need to be consistent with the row names of your
  data frames. It is good practice to use tip and node labels (or node
  numbers) when you merge data with a tree.

  If you provide both \code{tip.data} and \code{node.data}, the
  treatment of columns with common names will depend on the
  \code{merge.data} argument. If TRUE, columns with the same name in
  both data frames will be merged; when merging columns of different
  data types, coercion to a common type will follow standard R rules. If
  \code{merge.data} is FALSE, columns with common names will be
  preserved independently, with ".tip" and ".node" appended to the
  names. This argument has no effect if \code{tip.data} and
  \code{node.data} have no column names in common.

  If you provide \code{all.data} along with either of \code{tip.data}
  and \code{node.data}, it must have distinct column names, otherwise an
  error will result. Additionally, although supplying columns with the
  same names \emph{within} data frames is not illegal, automatic
  renaming for uniqeness may lead to surprising results, so this
  practice should be avoided.
}

\value{
  An object of class \linkS4class{phylo4d}.
}

\seealso{
\code{\link{coerce-methods}} for translation functions. The
  \linkS4class{phylo4d} class, the \code{\link{formatData}}
  function to check the validity of \code{phylo4d} objects;
  \linkS4class{phylo4} class and \link{phylo4} constructor.}

\author{Ben Bolker, Thibaut Jombart, Steve Kembel, Francois Michonneau,
  Jim Regetz}

\note{
  Checking on matches between the tree and the data will be done by the
  validity checker (label matches between data and tree tips, number of
  rows of data vs. number of nodes/tips/etc.)
}
\examples{
treeOwls <- "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);"
tree.owls.bis <- read.tree(text=treeOwls)
try(phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3)), silent=TRUE)
obj <- phylo4d(as(tree.owls.bis,"phylo4"),data.frame(wing=1:3), match.data=FALSE)
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
rNodeData <- data.frame(randomTrait = rnorm(nNodes(geoTree)),
                        row.names = nodeId(geoTree, "internal"))

exGeo2 <- phylo4d(geoTree, tip.data = geoData, node.data = rNodeData)

### Example using 'merge.data'
data(geospiza)
trGeo <- extractTree(geospiza)
tDt <- data.frame(a=rnorm(nTips(trGeo)), row.names=nodeId(trGeo, "tip"))
nDt <- data.frame(a=rnorm(nNodes(trGeo)), row.names=nodeId(trGeo, "internal"))

(matchData1 <- phylo4d(trGeo, tip.data=tDt, node.data=nDt, merge.data=FALSE))
(matchData2 <- phylo4d(trGeo, tip.data=tDt, node.data=nDt, merge.data=TRUE))

## Example with 'all.data'
nodeLabels(geoTree) <- as.character(nodeId(geoTree, "internal"))
rAllData <- data.frame(randomTrait = rnorm(nTips(geoTree) + nNodes(geoTree)),
row.names = labels(geoTree, 'all'))

exGeo5 <- phylo4d(geoTree, all.data = rAllData)

## Examples using 'rownamesAsLabels' and comparing with match.data=FALSE
tDt <- data.frame(x=letters[1:nTips(trGeo)],
                  row.names=sample(nodeId(trGeo, "tip")))
tipLabels(trGeo) <- as.character(sample(1:nTips(trGeo)))
(exGeo6 <- phylo4d(trGeo, tip.data=tDt, rownamesAsLabels=TRUE))
(exGeo7 <- phylo4d(trGeo, tip.data=tDt, rownamesAsLabels=FALSE))
(exGeo8 <- phylo4d(trGeo, tip.data=tDt, rownamesAsLabels=FALSE, match.data=FALSE))

require(ape) ## for rcoal
## generate a tree and some data
set.seed(1)
p3 <- rcoal(5)
dat <- data.frame(a = rnorm(5), b = rnorm(5), row.names = p3$tip.label)
dat.defaultnames <- dat
row.names(dat.defaultnames) <- NULL
dat.superset <- rbind(dat, rnorm(2))
dat.subset <- dat[-1, ]

## create a phylo4 object from a phylo object
p4 <- as(p3, "phylo4")

## create phylo4d objects with tip data
p4d <- phylo4d(p4, dat)
###checkData(p4d)
p4d.sorted <- phylo4d(p4, dat[5:1, ])
try(p4d.nonames <- phylo4d(p4, dat.defaultnames))
p4d.nonames <- phylo4d(p4, dat.defaultnames, match.data=FALSE)

\dontrun{
p4d.subset <- phylo4d(p4, dat.subset)
p4d.subset <- phylo4d(p4, dat.subset)
try(p4d.superset <- phylo4d(p4, dat.superset))
p4d.superset <- phylo4d(p4, dat.superset)
}

## create phylo4d objects with node data
nod.dat <- data.frame(a = rnorm(4), b = rnorm(4))
p4d.nod <- phylo4d(p4, node.data = nod.dat, match.data=FALSE)


## create phylo4 objects with node and tip data
p4d.all1 <- phylo4d(p4, node.data = nod.dat, tip.data = dat, match.data=FALSE)
nodeLabels(p4) <- as.character(nodeId(p4, "internal"))
p4d.all2 <- phylo4d(p4, all.data = rbind(dat, nod.dat, match.data=FALSE))


}
\keyword{misc}
