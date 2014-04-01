



setMethod("depthTips", signature(x="phylo4"), function(x) {
  nodeDepth(x, 1:nTips(x))
})



setMethod("nodeDepth", signature(x="phylo4"),
  function(x, node) {
    if (!hasEdgeLength(x))
      return(NULL)
    else {
      node <- getNode(x, node, missing="warn")
      node <- node[!is.na(node)]
      res <- sapply(node, function(n)
                    sumEdgeLength(x, ancestors(x, n, "ALL")))
      if (length(res) == 1) {
        res <- res[[1]]
        names(res) <- names(node)
      }      
      res
    }
})


setMethod("hasEdgeLength", signature(x="phylo4"),
 function(x) {
    !all(is.na(x@edge.length))
})

# return edge lengths in order by edgeIds (same order as edge matrix)
setMethod("edgeLength", signature(x="phylo4"),
 function(x, node) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    if (!missing(node)) {
        id <- getEdge(x, node)
    } else {
        id <- edgeId(x, "all")
    }
    elen <- x@edge.length[match(id, names(x@edge.length))]
    names(elen) <- id
    return(elen)
})

setReplaceMethod("edgeLength", signature(x="phylo4"),
 function(x, use.names=TRUE, ..., value) {
    len <- .createEdge(value, x@edge, type="lengths", use.names)
    ## return empty vector if all values are NA
    if (all(is.na(len))) {
        emptyVec <- numeric()
        attributes(emptyVec) <- list(names=character(0))
        x@edge.length <- emptyVec
    } else {
        x@edge.length <- len
    }
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
})

setMethod("sumEdgeLength", signature(x="phylo4"),
 function(x, node) {
    if(!hasEdgeLength(x))
        NULL
    else {
        nd <- getNode(x, node)
        iEdges <- which(x@edge[,2] %in% nd)
        sumEdges <- sum(x@edge.length[iEdges],na.rm=TRUE)
        sumEdges
    }
})

setMethod("isUltrametric", signature(x="phylo4"),
  function(x, tol=.Machine$double.eps^.5) {
    if (!hasEdgeLength(x)) {
      stop("The tree has no edge lengths.")
    }
    if (identical(all.equal.numeric(var(depthTips(x)), 0, tolerance=tol), TRUE)) {
      TRUE
    }
    else FALSE
  })
