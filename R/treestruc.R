## not bothering to check for zero branch lengths:
##   consensus is that this isn't very important,
##  and that it's simple enough to do
##   any(edgeLength(x)==0) if necessary
hasPoly <- function(object) {
  if(!checkPhylo4(object)) stop("to be used with a phylo4 object")
  if (nEdges(object)==0) return(FALSE)
  degree <- tabulate(edges(object, drop.root=TRUE)[, 1])
  any(degree > 2)
}



#' Test trees for polytomies, inline nodes, or reticulation
#' 
#' checks to see whether trees have (structural) polytomies, inline nodes
#' (i.e., nodes with a single descendant), or reticulation (i.e., nodes with
#' more than one ancestor)
#' 
#' 
#' @aliases hasSingle hasPoly hasRetic
#' @param object an object inheriting from class \code{phylo4}
#' @return Logical value
#' @note Some algorithms are unhappy with structural polytomies (i.e., >2
#' descendants from a node), with single-descendant nodes, or with
#' reticulation; these functions check those properties.  We haven't bothered
#' to check for zero branch lengths: the consensus is that it doesn't come up
#' much, and that it's simple enough to test \code{any(edgeLength(x) == 0)} in
#' these cases.  (Single-descendant nodes are used e.g. in OUCH, or in other
#' cases to represent events occurring along a branch.)
#' @author Ben Bolker
#' @keywords misc
#' @examples
#' 
#' tree.owls.bis <- ape::read.tree(text = "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);")
#' owls4 <- as(tree.owls.bis, "phylo4")
#' hasPoly(owls4)
#' hasSingle(owls4)
#' 
hasSingle <- function(object) {
  if(!checkPhylo4(object)) stop("to be used with a phylo4 object")
  if (nEdges(object)==0) return(FALSE)
  degree <- tabulate(edges(object, drop.root=TRUE)[, 1])
  any(degree == 1)
}

hasRetic <- function(object) {
  if(!checkPhylo4(object)) stop("to be used with a phylo4 object")
  if (nEdges(object)==0) return(FALSE)
  ancest <- tabulate(edges(object)[, 2])
  any(ancest > 1)
}


### TO BE FINISHED - Thibaut

# Returns a vector of logical 
# TRUE = this edge goes from an internal node to another
#internEdges <- function(object){
#  if(!checkPhylo4(object)) stop("to be used with a phylo4 object")
#  x <- object
#  isTips <- (tabulate(x@edge[,1]) == 0)
#  tips <- x@edge[isTips, 1]
#  inter <- is.na(match(x@edge[,2],tips))
#  return(inter)
#}

# Returns a vector of logical 
# TRUE = this edge goes from an internal node to a tip
#terminEdges <- function(object){
#  return(!internEdges(object))
#}

#isPoly <- function(object, position=c("all", "terminal", "internal")){
#  if(!checkPhylo4(object)) stop("to be used with a phylo4 object")
#  x <- object
#  pos <- match.arg(position)
#  res <- (tabulate(x@edge[,1]) > 2)

  # all polytomies 
#  if(pos=="all") return(res)

  # find which edge ends at a tip
  
  
  
  # external polytomies
  
#}
