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
  \item{tip.data}{a data frame for tips data}
  \item{node.data}{a data frame for nodes data}
  \item{all.data}{a data frame for all (i.e. tips and nodes) data. In such
    case, first rows should correspond to tips, last rows to nodes.}
  \item{match.data}{(logical) should the rownames of the data frame provided in
  \code{tip.data}, \code{node.data} and/or \code{all.data} be used to be matched against
  tip and node labels?}
  \item{merge.data}{if both \code{tip.data} and \code{node.data} are
    provided, it determines if they should be merged as a single
    trait (if TRUE). This argument is evaluated only if both \code{tip.data} and
    \code{node.data} have identical column names.}
  \item{check.node.labels}{if \code{x} is of class \code{phylo}, use
    either "keep" (the default) to retain internal node labels, "drop"
    to drop them, or "asdata" to convert them to numeric tree data. This
    argument is useful if the \code{phylo} object has non-unique node
    labels or node labels with informative data (e.g., posterior
    probabilities).}
  \item{rownamesAsLabels}{(logical), in the case of number-like labels 
    should the row names of the data provided be considered as labels 
    (TRUE) or node numbers (FALSE and default)}
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
  tips and/or nodes. If you provide \code{all.data} and \code{tip.data}
  or \code{node.data}, row names of the data frames will be matched
  (\code{all.data} names are matched against \code{tip.data} and/or
  \code{node.data}). This is done independently of the labels of the tree
  (and also of the value of the arguments \code{use.tip.names} and
  \code{use.node.names}). This means that you need to be consistent
  with the row names of your data frames. It is good practice to use tip
  and node labels (or node numbers) when you merge data with a tree.
  If \code{all.data} is provided and tip and node names are not being
  used, then the first part of the data frame will be tips data and the
  second part will be nodes data.
}

\value{
  An object of class \linkS4class{phylo4d}.
}

\seealso{
\code{\link{coerce-methods}} for translation functions. The
  \linkS4class{phylo4d} class, the \code{\link{formatData}}
  function to check the validity of \code{phylo4d} objects;
  \linkS4class{phylo4} class and \link{phylo4} constructor.}

\author{Ben Bolker, Thibaut Jombart, Steve Kembel, Francois Michonneau}

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
