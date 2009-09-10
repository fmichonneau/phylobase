
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
