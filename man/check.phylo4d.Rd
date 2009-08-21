\name{checkPhylo4d}
\alias{formatData}

\alias{checkPhylo4Data}
\title{Data addition and validity checking for phylo4d objects}
\description{
  Basic checks on the validity of S4 phylogenetic plus data objects
  (i.e. \linkS4class{phylo4d} objects).
}
\section{Methods}{
  \describe{
	\item{formatData}{\code{signature(phy="phylo4d")}} 
	\item{checkPhylo4Data}{\code{signature(phy="phylo4d")}}
  }
  }

  \usage{
formatData(phy, dt, type=c("tip", "internal", "all"),
                       match.data=TRUE, rownamesAsLabels=FALSE,
		       label.type=c("rownames", "column"),
                       label.column=1, missing.data=c("fail", "warn", "OK"),
                       extra.data=c("warn", "OK", "fail"))
checkPhylo4Data(phy)
                       

}
\arguments{
  \item{phy}{A phylo4 object}
  \item{dt}{a data frame or matrix}
  \item{type}{type of data to attach}
  \item{match.data}{(logical) check data labels against tree labels?}
  \item{label.type}{character, "rownames" or "column": should
    the labels be taken from the row names of \code{dt} or from
    the \code{label.column} column of \code{dt}?}
  \item{label.column}{if \code{label.type=="column"}, column specifier
  (number or name) of the column containing tip labels}
  \item{missing.data}{action to take if there are missing data/labels don't match}
  \item{extra.data}{action to take if there are extra data/labels don't match} 
  \item{rownamesAsLabels}{(logical), in the case of number-like labels should the row
  names of the data provided be considered as labels (TRUE) or node numbers (FALSE and default)}
  
}
\value{
By default, tip data names are used to attach the tip data to the tree, must be present, and must match the set of tree taxa exactly. By default node data names are not used to attach the data to the tree, and thus node data names are not required.


}
\details{
  \code{formatData} is an internal function that should not be called by the user. It
  formats properly the data provided by the 
  user before being attached to a tree. The function (1) converts labels provided in the
  data into node numbers, (2) makes sure that the data are appropriately matched against
  tip or nodes, (3) checks for differences between data and tree, (4) creates a data frame
  with the correct dimensions given a tree. \code{formatData} returns a data frame. 
  
  \code{checkPhylo4Data} checks the validity of a \code{phylo4d} object. It ensures that
  (1) the data associated with the tree have the correct dimensions, (2) that the row
  names for the data are correct.
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
nodeLabels(p4) <- as.character(sort(nodeId(p4)))
p4d.all2 <- phylo4d(p4, all.data = rbind(dat, nod.dat, match.data=FALSE))

## examples of merge.data
data(geospiza)
trGeo <- extractTree(geospiza)
tDt <- data.frame(a=rnorm(nTips(trGeo)), row.names=nodeId(trGeo, "tip"))
nDt <- data.frame(a=rnorm(nNodes(trGeo)), row.names=nodeId(trGeo, "internal"))

(matchData1 <- phylo4d(trGeo, tip.data=tDt, node.data=nDt, merge.data=FALSE))
(matchData2 <- phylo4d(trGeo, tip.data=tDt, node.data=nDt, merge.data=TRUE))

}
\keyword{misc}
