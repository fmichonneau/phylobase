\name{checkPhylo4d}
\alias{checkData}
\alias{attachData}
\title{Validity checking for phylo4d objects}
\description{
  Basic checks on the validity of S4 phylogenetic plus data objects
  (i.e. \linkS4class{phylo4d} objects).
}
\usage{
checkData(object, label.type = c("row.names", "column"),
label.column = 1, use.tip.names=TRUE, missing.tip.data=c("fail","OK","warn"),
extra.tip.data=c("fail","OK","warn"), default.tip.names=c("warn","OK","fail"),
use.node.names=FALSE, missing.node.data=c("OK","warn","fail"),
extra.node.data=c("OK","warn","fail"), default.node.names=c("warn","OK","fail"),
non.unique.tips=c("warn", "OK", "fail"), non.unique.nodes=c("warn", "OK", "fail"), \dots)

attachData(object, label.type = c("row.names", "column"), label.column = 1,
use.tip.names = TRUE, use.node.names = FALSE,\dots)
}

\arguments{
  \item{object}{A phylo4d object}
  \item{label.type}{Obtain labels for matching data to tree labels from the row.names or a column of the data? (default=\code{row.names})}
  \item{label.column}{If label.type="column", number or name of column to use for matching data to tree labels (default=\code{1})}
  \item{use.tip.names}{Use tip data names if present (default=\code{TRUE})}
  \item{missing.tip.data}{Can tip data taxa be a subset of tree taxa? (default=\code{fail})}
  \item{extra.tip.data}{Can tip data taxa be a superset of tree taxa? (default=\code{fail})}
  \item{default.tip.names}{Do tip data have the default (uninformative) names? (default=\code{warn})}
  \item{use.node.names}{Use node data names if present
    (default=\code{FALSE})}
  \item{missing.node.data}{Can node data taxa be a subset of tree taxa?
    (default=\code{OK})}
  \item{extra.node.data}{Can node data taxa be a superset of tree taxa?
    (default=\code{OK})}
  \item{default.node.names}{Do tip data have the default (uninformative) names?
  (default=\code{warn})}
  \item{non.unique.tips}{Can you attach data to non-unique tip? (default=\code{warn})}
  \item{non.unique.nodes}{Can you attach data to non-unique nodes? (default=\code{warn})}
  \item{\dots}{Further arguments (for consistency)}
}
\value{
  \code{checkData} returns an error string (describing problems) or \code{TRUE} if everything is OK. If option \code{fail} is specified and a problem is encountered, returns an error string. If option \code{warn} is specified and a problem is encountered, returns \code{TRUE} plus a warning message.

  By default, tip data names are used to attach the tip data to the tree, must be present, and must match the set of tree taxa exactly. By default node data names are not used to attach the data to the tree, and thus node data names are not required.

  \code{attachData} returns the phylo4 object with tipdata and nodedata (if present) attached to the phylo4 tree.
}
\details{
  \code{checkData} checks validity of a phylo4d object. For tipdata and nodedata, checks
  for presence of tip and node names, whether tip or node data taxa and tree taxa are an
  exact match, and whether tip or node data taxa are a superset or subset of the tree
  taxa. By default, if several tip or node labels are non-unique, they all receive the same
  value. This is only true if you use names to attach data to your tree. You can prevent
  this behavior by changing the value of the arguments \code{non.unique.tips} and
  \code{non.unique.nodes} to \code{fail}.

  \code{attachData} takes a phylo4d object that has been checked, and attaches the tipdata and nodedata to the object. This involves the data being subset to include only the taxa present in the tree, with missing data taxa assigned NA values. The data are sorted to match the order of the tip.labels of the tree to allow use with comparative analysis functions.
}
\seealso{the \code{\link{phylo4d}} constructor, the \linkS4class{phylo4d} class. See also the \code{\link{checkPhylo4}}, the \code{\link{phylo4}} constructor and the \linkS4class{phylo4} class. See \code{\link{coerce-methods}} for translation functions.}
\author{Steven Kembel \email{skembel@berkeley.edu}}
\examples{
require(ape)
#generate a tree and some data
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
