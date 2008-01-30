\name{extract.tree}
\alias{extract.tree}
\title{Get tree from tree+data object}
\description{
  Extracts a \code{phylo4} tree object from a \code{phylo4d} tree+data object.
}
\usage{
extract.tree(from)
}
\arguments{
  \item{from}{ a \code{phylo4d} object, containing a phylogenetic tree plus associated phenotypic data. Created from \code{phylo4d() function}. }
}
\details{
  \code{extract.tree} extracts just the phylogeny from a tree+data object. The phylogeny contains the topology (how the nodes are linked together), the branch lengths (if any), and any tip and/or node labels. This may be useful for extracting a tree from a \code{phylo4d} object, and associating with another phenotypic dataset, or to convert the tree to another format. 
}
\author{ Ben Bolker }
\seealso{\code{\link{phylo4}}, \code{\link{phylo4d}}, \code{\link{coerce-methods}} for translation functions. }
\examples{
tree.phylo <- read.tree(text="((a,b),c);")
tree <- as(tree.phylo, "phylo4")
plot(tree)
tip.data <- data.frame(size=c(1,2,3), row.names=c("a", "b", "c"))
treedata <- phylo4d(tree, tip.data)
plot(treedata)
print(treedata)

tree1 <- extract.tree(treedata)
print(tree1)
}
\keyword{methods}
