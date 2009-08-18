\name{checkPhylo4d}
\alias{formatData}
\alias{addData,phylo4d-method}
\title{Data addition and validity checking for phylo4d objects}
\description{
  Basic checks on the validity of S4 phylogenetic plus data objects
  (i.e. \linkS4class{phylo4d} objects).
}
\section{Methods}{
  \describe{
    \item{addData}{\code{signature(x = "phylo4d")}: add data}
  }
  }

  \usage{
formatData(phy, dt, type=c("tip", "internal", "all"),
                       match.data=TRUE, label.type=c("rownames", "column"),
                       label.column=1, missing.data=c("fail", "warn", "OK"),
                       extra.data=c("warn", "OK", "fail"))
\S4method{addData}{phylo4d}(x, tip.data=NULL, node.data=NULL,
                                         all.data=NULL, pos=c("after", "before"),
                                         merge.data=TRUE, match.data=TRUE,
                                         \dots)
}
\arguments{
  \item{x}{A phylo4d object}
  \item{phy}{A phylo4 object}
  \item{dt}{a data frame or matrix}
  \item{type}{type of data to attach}
  \item{label.type}{stuff}
  \item{label.column}{stuff}
  \item{missing.data}{action to take if there are missing data/labels don't match}
  \item{extra.data}{action to take if there are extra data/labels don't match}
  \item{tip.data}{a data frame, matching tree tips}
  \item{node.data}{a data frame, matching tree internal nodes}
  \item{all.data}{a data frame, matching all tree nodes (tips and internal nodes)}
  \item{pos}{add data columns before or after existing data?}
  \item{merge.data}{stuff}
  \item{match.data}{stuff}
  \item{\dots}{stuff}
}
\value{
  By default, tip data names are used to attach the tip data to the tree, must be present, and must match the set of tree taxa exactly. By default node data names are not used to attach the data to the tree, and thus node data names are not required.

  \code{addData} returns the phylo4 object with tipdata and nodedata (if present) attached to the phylo4 tree.
}
\details{
  \code{formatData} checks validity of a phylo4d object. For tipdata and nodedata, checks
  for presence of tip and node names, whether tip or node data taxa and tree taxa are an
  exact match, and whether tip or node data taxa are a superset or subset of the tree
  taxa. By default, if several tip or node labels are non-unique, they all receive the same
  value. This is only true if you use names to attach data to your tree. You can prevent
  this behavior by changing the value of the arguments \code{non.unique.tips} and
  \code{non.unique.nodes} to \code{fail}.
}
\seealso{the \code{\link{phylo4d}} constructor, the \linkS4class{phylo4d} class. See also the \code{\link{checkPhylo4}}, the \code{\link{phylo4}} constructor and the \linkS4class{phylo4} class. See \code{\link{coerce-methods}} for translation functions.}
\author{Steven Kembel \email{skembel@berkeley.edu}}
\examples{
require(ape) ## for rcoal
## generate a tree and some data
set.seed(1)
p3 <- rcoal(5)
dat <- data.frame(a = rnorm(5), b = rnorm(5), row.names = p3$tip.label)
dat.defaultnames <- dat
row.names(dat.defaultnames) <- NULL
dat.superset <- rbind(dat, rnorm(2))
dat.subset <- dat[-1, ]

#create a phylo4 object from a phylo object
p4 <- as(p3, "phylo4")

#create phylo4d objects with tip data
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

#create phylo4d objects with node data
nod.dat <- data.frame(a = rnorm(4), b = rnorm(4))
p4d.nod <- phylo4d(p4, node.data = nod.dat, match.data=FALSE)


#create phylo4 objects with node and tip data
p4d.all1 <- phylo4d(p4, node.data = nod.dat, tip.data = dat, match.data=FALSE)
nodeLabels(p4) <- as.character(sort(nodeId(p4)))
p4d.all2 <- phylo4d(p4, all.data = rbind(dat, nod.dat, match.data=FALSE))
}
\keyword{misc}
