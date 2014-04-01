
##' Number of tips, nodes and edges found in a tree.
##'
##' Function to return the number of tips, nodes and edges found in a
##' tree in the \code{phylo4} or \code{phylo4d} format.
##' @title nTips, nNodes, nEdges
##' @param x a \code{phylo4} or \code{phylo4d} object
##' @return a numeric vector indicating the number of tips, nodes or
##' edge respectively.
##' @docType methods
##' @export
##' @include phylo4-methods.R
##' @rdname nTips-methods
setGeneric("nTips", function(x) {
    standardGeneric("nTips")
})

##' @rdname nTips-methods
##' @aliases nTips,phylo4-method
setMethod("nTips", signature(x="phylo4"), function(x) {
    E <- edges(x)
    if(nrow(E) == 0)
        return(0)
    else {
        ## at this time NAs are not allowed in edge matrix
        ## sum(tabulate(E[, 1]) == 0)
        nTipsFastCpp(E[, 1])
    }        
})

##' @rdname nTips-methods
##' @aliases nTips,phylo-method
setMethod("nTips", signature(x="phylo"),
 function(x) {
     Ntip(x)
})

##' @rdname nTips-methods
##' @aliases nNodes
setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})

##' @rdname nTips-methods
##' @aliases nNodes,phylo4-method
setMethod("nNodes", signature(x="phylo4"), function(x) {
    E <- edges(x, drop.root=TRUE)
    if(nrow(E) == 0) {
        return(0)
    } else {
        return(length(unique(E[, 1])))
    }
})

##' @rdname nTips-methods
##' @aliases nEdges
setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})

##' @rdname nTips-methods
##' @aliases nEdges,phylo4-method
setMethod("nEdges", signature(x="phylo4"),
 function(x) {
    nrow(x@edge)
})


#########################################################
### Edge accessors
#########################################################

##' edges accessors
##'
##' @param x A \code{phylo4} or \code{phylo4d} object.
##'
##' @return \code{edges} returns the edge matrix that represent the
##' ancestor-descendant relationships among the nodes of the tree.
##'
##' \code{edgeOrder} returns the order in which the edge matrix is in.
##'
##' @seealso reorder
##' @export
##' @docType methods
##' @rdname edges-accessors
##' @include phylo4-methods.R
setGeneric("edges", function(x, ...) {
    standardGeneric("edges")
})

##' @rdname edges-accessors
##' @aliases edges,phylo4-method
setMethod("edges", signature(x="phylo4"),
 function(x, drop.root=FALSE, ...) {
     e <- x@edge
     if (drop.root) e <- e[e[, 1] != 0, ]
     e
})

##' @rdname edges-accessors
##' @aliases edgeOrder,phylo4-method
setMethod("edgeOrder", signature(x="phylo4"),
 function(x, ...) {
    x@order
})

