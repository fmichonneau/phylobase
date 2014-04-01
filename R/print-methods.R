
#' print a phylogeny
#' 
#' Prints a phylo4 or phylo4d object in data.frame format with user-friendly
#' column names
#' 
#' This is a user-friendly version of the tree representation, useful for
#' checking that objects were read in completely and translated correctly. The
#' phylogenetic tree is represented as a list of numbered nodes, linked in a
#' particular way through time (or rates of evolutionary change).  The topology
#' is given by the pattern of links from each node to its ancestor. Also given
#' are the taxon names, node type (root/internal/tip) and phenotypic data (if
#' any) associated with the node, and the branch length from the node to its
#' ancestor. A list of nodes (descendants) and ancestors is minimally required
#' for a phylo4 object.
#' 
#' @param x a \code{phylo4} tree or \code{phylo4d} tree+data object
#' @param edgeOrder in the data frame returned, the option 'pretty' returns the
#' internal nodes followed by the tips, the option 'real' returns the nodes in
#' the order they are stored in the edge matrix.
#' @param printall default prints entire tree. printall=FALSE returns the first
#' 6 rows
#' @return A data.frame with a row for each node (descendant), sorted as
#' follows: root first, then other internal nodes, and finally tips.\cr The
#' returned data.frame has the following columns:\cr \item{label}{Label for the
#' taxon at the node (usually species name).} \item{node}{Node number, i.e. the
#' number identifying the node in edge matrix.} \item{ancestor}{Node number
#' of the node's ancestor.} \item{branch.length}{The branch length connecting
#' the node to its ancestor (NAs if missing).} \item{node.type}{"root",
#' "internal", or "tip". (internally generated)} \item{data}{phenotypic data
#' associated with the nodes, with separate columns for each variable.}
#' @note This is the default show() method for phylo4, phylo4d. It prints the
#' user-supplied information for building a phylo4 object. For a full
#' description of the phylo4 S4 object and slots, see \code{\link{phylo4}}.
#' @author Marguerite Butler Thibaut Jombart
#' \email{jombart@@biomserv.univ-lyon1.fr} Steve Kembel
#' @keywords methods
#' @examples
#' 
#' 
#' tree.phylo <- ape::read.tree(text="((a,b),c);")
#' tree <- as(tree.phylo, "phylo4")
#' ##plot(tree,show.node=TRUE) ## plotting broken with empty node labels: FIXME
#' tip.data <- data.frame(size=c(1,2,3), row.names=c("a", "b", "c"))
#' treedata <- phylo4d(tree, tip.data)
#' plot(treedata)
#' print(treedata)
#' 
#' 
#' @export printphylo4
printphylo4 <- function(x, edgeOrder=c("pretty", "real"), printall=TRUE) {
    if(!nrow(edges(x))) {
        msg <- paste("Empty \'", class(x), "\' object\n", sep="")
        cat(msg)
    }
    else {
        toRet <- .phylo4ToDataFrame(x, edgeOrder)
        if (printall) {
            print(toRet)
        }
        else {
            print(head(toRet))
        }
    }
}

### Hack for print/show
### from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html
setMethod("print", "phylo4", printphylo4)
setMethod("show", signature(object="phylo4"),
   function(object) printphylo4(object))

### names
setMethod("names", signature(x="phylo4"),
 function(x) {
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})

### Head and Tail
setMethod("head", signature(x="phylo4"),
  function(x, n=20) {
      head(as(x,"data.frame"),n=n)
  })

setMethod("tail", signature(x="phylo4"),
  function(x, n=20) {
      tail(as(x, "data.frame"), n=n)
  })

