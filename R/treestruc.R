
## not bothering to check for zero branch lengths:
##   consensus is that this isn't very important,
##  and that it's simple enough to do
##   any(edgeLength(x)==0) if necessary
hasPoly <- function(object) {
  if(!check_phylo4(object)) stop("to be used with a phylo4 object")
  degree <- tabulate(edges(object)[, 1])
  struc <- any(degree > 2)
  return(struc)
}



hasSingles <- function(object) {
  if(!check_phylo4(object)) stop("to be used with a phylo4 object")
  degree <- tabulate(edges(object)[, 1])
  any(degree == 1)
}

hasRetic <- function(object) {
  if(!check_phylo4(object)) stop("to be used with a phylo4 object")
  ancest <- tabulate(edges(object)[, 2])
  any(ancest>1)
}


### TO BE FINISHED - Thibaut

# Returns a vector of logical 
# TRUE = this edge goes from an internal node to another
#internEdge <- function(object){
#  if(!check_phylo4(object)) stop("to be used with a phylo4 object")
#  x <- object
#  isTips <- (tabulate(x@edge[,1]) == 0)
#  tips <- x@edge[isTips, 1]
#  inter <- is.na(match(x@edge[,2],tips))
#  return(inter)
#}

# Returns a vector of logical 
# TRUE = this edge goes from an internal node to a tip
#terminEdge <- function(object){
#  return(!internEdge(object))
#}

#isPoly <- function(object, position=c("all", "terminal", "internal")){
#  if(!check_phylo4(object)) stop("to be used with a phylo4 object")
#  x <- object
#  pos <- match.arg(position)
#  res <- (tabulate(x@edge[,1]) > 2)

  # all polytomies 
#  if(pos=="all") return(res)

  # find which edge ends at a tip
  
  
  
  # external polytomies
  
#}
