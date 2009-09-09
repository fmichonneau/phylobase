### This file contains the methods and accessors for phylo4(d) objects
### The file is organized in sections:

### 1. Tip accessors
###  1.1. nTips()

### 2. Node accessors
###  2.1. nNodes()
###  2.2. nodeType()
###  2.3. nodeId()

### 3. Edge accessors
###  3.1. nEdges()
###  3.2. edges()
###  3.3. edgeOrder()
###  3.4. edgeId()
###  3.5. hasEdgeLength()
###  3.6. edgeLength()
###  3.7. edgeLength() <-
###  3.8. sumEdgeLength()

### 4. Root accessors
###  4.1. isRooted()
###  4.2. rootNode()
###  4.3. rootNode() <-

### 5. Label accessors
###  5.1. labels()
###  5.2. labels() <-
###  5.3. hasNodeLabels()
###  5.4. nodeLabels()
###  5.5. nodeLabels() <-
###  5.6. tipLabels()
###  5.7. tipLabels() <-
###  5.8. hasEdgeLabels()
###  5.9. edgeLabels()
###  5.10. edgeLabels() <-

### 6. Displaying functions
###  6.1. printphylo4()
###  6.2. print()
###  6.3  show()
###  6.4. names()
###  6.5. head()
###  6.6. tail()
###  6.7. summary()

### 7. Ordering
###  7.1. orderIndex()
###  7.2. reorder()


#########################################################
### Tip accessors
#########################################################

setMethod("nTips", "phylo4", function(x, ...) {
    E <- edges(x)
    if(nrow(E) == 0)
        return(0)
    else {
        ## doesn't handle reticulated networks
        ##    res <- sum(!E[, 2] %in% E[, 1])
        res <- sum(tabulate(na.omit(E[,1])) == 0) ## twice as fast as ...
        ## change suggested by Aaron Mackey, handles reticulated networks better
        ## res <- sum(!(unique(E[,2]) %in% E[,1]))
        return(res)
    }
})

## hack to ensure ape compatibility
setMethod("nTips","ANY", function(x) {
    if (class(x)=="phylo") {
        Ntip(x)
    } else stop(paste("no 'nTips' method available for",
                      deparse(substitute(x)),
                      "(class",class(x),")"))
})

#########################################################
### Node accessors
#########################################################

setMethod("nNodes", "phylo4", function(x) {
    x@Nnode
})

setMethod("nodeType", "phylo4", function(phy) {
    if(nTips(phy) == 0)
        return(NULL)
    else {
        listNodes <- sort(unique(as.vector(edges(phy))))
        t <- rep("internal", length(listNodes)) # FM: internal is default (I think it's safer)
        names(t) <- listNodes

        ## node number of real internal nodes
        iN <- names(table(edges(phy)[,1]))
        ## node number that are not internal nodes (ie that are tips)
        tN <- names(t)[!names(t) %in% iN]
        t[tN] <- "tip"

        ## if the tree is rooted
        if(isRooted(phy)) t[rootNode(phy)] <- "root"

        return(t)
    }
})

# return node IDs (or a subset thereof) in ascending order
setMethod("nodeId", "phylo4", function(x, type=c("all",
    "tip","internal","root")) {

     type <- match.arg(type)
     E <- edges(x)

     ## Note: this implementation will still work even if tips are not
     ## 1:nTips and nodes are not (nTips+1):nNodes
     nid <- switch(type,
         ## all nodes appear at least once in the edge matrix
         all = unique(na.omit(as.vector(E))),
         ## tips are nodes that do not appear in the ancestor column
         tip = setdiff(E[, 2], E[, 1]),
         ## internals are nodes that *do* appear in the ancestor column
         internal = na.omit(unique(E[, 1])),
         ## roots are nodes that have NA as ancestor
         root = if (!isRooted(x)) NA else unname(E[is.na(E[, 1]), 2]))

     return(sort(nid))

})



#########################################################
### Edge accessors
#########################################################

setMethod("nEdges", "phylo4", function(x) {
    nrow(x@edge)
})

# return edge matrix in its current order
setMethod("edges", "phylo4", function(x, order, drop.root=FALSE, ...) {
  e <- x@edge
  if (drop.root) e <- e[!is.na(e[,1]),]
  e
})

setMethod("edgeOrder", "phylo4", function(x, ...) {
    x@order
})

# return edge IDs (or a subset thereof) in edge matrix order
setMethod("edgeId", "phylo4", function(x, type=c("all", "tip",
    "internal", "root")) {
    type <- match.arg(type)
    edge <- edges(x)
    if (type=="tip") {
        isTip <- !(edge[, 2] %in% edge[, 1])
        edge <- edge[isTip, , drop=FALSE]
    } else if (type=="internal") {
        isInt <- (edge[, 2] %in% edge[, 1])
        edge <- edge[isInt, , drop=FALSE]
    } else if (type=="root") {
        isRoot <- is.na(edge[, 1])
        edge <- edge[isRoot, , drop=FALSE]
    } # else just use complete edge matrix if type is "all"
    id <- paste(edge[, 1], edge[, 2], sep="-")
    return(id)
})

setMethod("hasEdgeLength","phylo4", function(x) {
    !all(is.na(x@edge.length))
})

# return edge lengths in order by edgeIds (same order as edge matrix)
setMethod("edgeLength", "phylo4", function(x, node) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    elen <- x@edge.length[match(edgeId(x, "all"), names(x@edge.length))]
    if (!missing(node)) {
        n <- getNode(x, node)
        elen <- elen[match(n, x@edge[,2])]
    }
    return(elen)
})

setReplaceMethod("edgeLength", "phylo4", function(x, use.names=TRUE, ..., value) {
    if(use.names && !is.null(names(value))) {
        if(!all(names(value) %in% names(x@edge.length)))
            stop("Names provided don't match internal edge labels")
        x@edge.length[match(names(value), names(x@edge.length))] <- value
    }
    else
        x@edge.length[1:nEdges(x)] <- value
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
})

setMethod("sumEdgeLength", "phylo4", function(phy, node) {
    if(!hasEdgeLength(phy))
        NULL
    else {
        nd <- getNode(phy, node)
        iEdges <- which(phy@edge[,2] %in% nd)
        sumEdges <- sum(phy@edge.length[iEdges],na.rm=TRUE)
        sumEdges
    }
})

#########################################################
### Root accessors
#########################################################

setMethod("isRooted","phylo4", function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x) == 0) return(FALSE)
    any(is.na(edges(x)[,1]))
})

setMethod("rootNode", "phylo4", function(x) {
    if (!isRooted(x))
        return(NA)
    unname(edges(x)[which(is.na(edges(x)[,1])),2])
})

setReplaceMethod("rootNode", "phylo4", function(x, value) {
    stop("Root node replacement not implemented yet")
})

#########################################################
### Label accessors
#########################################################

## return labels in increasing node order
setMethod("labels", "phylo4", function(object, type = c("all", "tip",
    "internal")) {
    type <- match.arg(type)
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    if (type=="all") {
        all <- c(object@tip.label, object@node.label)
        return(all[match(nodeId(object, "all"), names(all))])
    } else if (type=="tip") {
        tip <- object@tip.label
        return(tip[match(nodeId(object, "tip"), names(tip))])
    } else if (type=="internal") {
        int <- object@node.label
        return(int[match(nodeId(object, "internal"), names(int))])
    }
})

setReplaceMethod("labels",
                 signature(object="phylo4", type="ANY",
                           use.names="ANY", value="character"),
   function(object, type = c("tip", "internal", "allnode"),
            use.names, ..., value) {

       ## Default options
       if(missing(type))
           type <- "tip"
       if (missing(use.names))
           use.names <- FALSE

       type <- match.arg(type)


       ob <- switch(type,
              ## If 'tip'
              tip = {
                  object@tip.label <- .createLabels(value, nTips(object),
                                                    nNodes(object), use.names,
                                                    type="tip")
                  object
              },
              ## If 'internal'
              internal = {
                  object@node.label <- .createLabels(value, nTips(object),
                                                     nNodes(object), use.names,
                                                     type="internal")
                  object
              },
              ## If 'allnode'
              allnode = {
                  if(use.names) {
                      tipVal <- value[names(value) %in% nodeId(object, "tip")]
                      nodVal <- value[names(value) %in% nodeId(object, "internal")]
                      object@tip.label <- .createLabels(tipVal, nTips(object),
                                                        nNodes(object), use.names,
                                                        type="tip")
                      object@node.label <- .createLabels(nodVal, nTips(object),
                                                         nNodes(object), use.names,
                                                         type="internal")
                  }
                  else {
                      ntips <- nTips(object)
                      nedges <- nTips(object) + nNodes(object)
                      object@tip.label <- .createLabels(value[1:ntips], nTips(object),
                                                        nNodes(object), use.names,
                                                        type="tip")
                      object@node.label <- .createLabels(value[(ntips+1):nedges],
                                                         nTips(object),
                                                         nNodes(object), use.names,
                                                         type="internal")
                  }
                  object
              })

       if(is.character(checkval <- checkPhylo4(ob)))
           stop(checkval)
       else
           return(ob)
   })


### Node Labels
setMethod("hasNodeLabels", "phylo4", function(x) {
    !all(is.na(x@node.label))
})

setMethod("nodeLabels", "phylo4", function(object) {
    labels(object, type="internal")
})

setReplaceMethod("nodeLabels", signature(object="phylo4", value="character"),
  function(object, ..., value) {
      labels(object, type="internal", ...) <- value
      if(is.character(checkval <- checkPhylo4(object))) stop(checkval)
      object
  })

### Tip labels
setMethod("tipLabels", "phylo4", function(object) {
    labels(object, type="tip")
    })

setReplaceMethod("tipLabels", signature(object="phylo4", value="character"),
  function(object, ...,  value) {
      labels(object, type="tip", ...) <- value
      if(is.character(checkval <- checkPhylo4(object))) stop(checkval)
      return(object)
  })


### Edge labels
setMethod("hasEdgeLabels", "phylo4", function(x) {
    length(x@edge.label) > 0
})

# return edge labels in order by edgeIds (same order as edge matrix)
setMethod("edgeLabels", signature(x = "phylo4"), function(x) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    x@edge.label[match(edgeId(x, "all"), names(x@edge.label))]
})

setReplaceMethod("edgeLabels", signature(object="phylo4", value="character"),
  function(object, ..., value) {
      object@edge.label <- value
      if(is.character(checkval <- checkPhylo4(object))) stop(checkval)
      object
  })


#########################################################
### Displaying phylo4: print, show, head, tail and summary
#########################################################

### print
printphylo4 <- function(x, edgeOrder=c("pretty", "real"), printall=TRUE) {
    if(!nrow(x@edge)) {
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
setMethod("show", "phylo4", function(object) printphylo4(object))

### names
setMethod("names", signature(x = "phylo4"), function(x){
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})

### Head and Tail
setMethod("head",signature(x = 'phylo4'),
          function(x,n=20) {
            head(as(x,"data.frame"),n=n)
          })

setMethod("tail",signature(x = 'phylo4'),
          function(x,n=20) {
            tail(as(x,"data.frame"),n=n)
          })

### summary
setMethod("summary","phylo4", function (object, quiet=FALSE) {
    x <- object
    res <- list()

    ## build the result object
    res$name <- deparse(substitute(object, sys.frame(-1)))
    res$nb.tips <- length(x@tip.label)
    res$nb.nodes <- x@Nnode

    if(!is.null(x@edge.length)){
        res$mean.el <- mean(x@edge.length, na.rm=TRUE)
        res$var.el <- var(x@edge.length, na.rm=TRUE)
        if (isRooted(x) && is.na(x@edge.length[rootNode(x)])) {
            res$sumry.el <- summary(x@edge.length[-rootNode(x)])
        } else {
            res$sumry.el <- summary(x@edge.length)
        }
    } else {
        res$mean.el <- NULL
        res$var.el <- NULL
        res$sumry.el <- NULL
    }

    ## check for polytomies
    if (nrow(edges(x)) != 0 && any(tabulate(na.omit(edges(object)[,1]))>2)){ # if there are polytomies
        E <- edges(x)
        temp <- tabulate(na.omit(E[,1]))
        degree <- temp[na.omit(E[,1])] # contains the degree of the ancestor for all edges
        endsAtATip <- !(E[,2] %in% E[,1])
        terminPoly <- (degree>2) & endsAtATip
        internPoly <- (degree>2) & !endsAtATip
        res$degree <- degree
        res$polytomy <- rep("none",nrow(E))
        res$polytomy[terminPoly] <- "terminal"
        res$polytomy[internPoly] <- "internal"
        ## now just keep information about nodes (not all edges)
        nod <- unique(E[,1])
        idx <- match(nod,E[,1])
        res$degree <- res$degree[idx]
        names(res$degree) <- nodeLabels(x)
        res$polytomy <- res$polytomy[idx]
        names(res$polytomy) <- nodeLabels(x)
    }

    ## model info
    res$loglik <- attr(x, "loglik")
    res$para <- attr(x, "para")
    res$xi <- attr(x, "xi")

    ## if quiet, stop here
    if(quiet) return(invisible(res))

    ## now, print to screen is !quiet
    cat("\n Phylogenetic tree :", res$name, "\n\n")
    cat(" Number of tips    :", res$nb.tips, "\n")
    cat(" Number of nodes   :", res$nb.nodes, "\n")
    ## cat("  ")
    if(!length(x@edge.length)) {
        cat(" Branch lengths    : No branch lengths.\n")
    } else {
        cat(" Branch lengths:\n")
        cat("        mean         :", res$mean.el, "\n")
        cat("        variance     :", res$var.el, "\n")
        cat("        distribution :\n")
        print(res$sumry.el)
    }
    if(nrow(edges(x)) != 0 && hasPoly(x)){
        cat("\nDegree of the nodes  :\n")
        print(res$degree)
        cat("\n")
        cat("Types of polytomy:\n")
        print(res$polytomy)
        cat("\n")
    }

    if (!is.null(attr(x, "loglik"))) {
        cat("Phylogeny estimated by maximum likelihood.\n")
        cat("  log-likelihood:", attr(x, "loglik"), "\n\n")
        npart <- length(attr(x, "para"))
        for (i in 1:npart) {
            cat("partition ", i, ":\n", sep = "")
            print(attr(x, "para")[[i]])
            if (i == 1)
              next
            else cat("  contrast parameter (xi):", attr(x,"xi")[i - 1], "\n")
        }
    }
    return(invisible(res))

}) # end setMethod summary phylo4




#########################################################
### Ordering
#########################################################

orderIndex <- function(phy, order = c('preorder', 'postorder')) {

    order <- match.arg(order)

    ## get a root node free edge matrix
    edge <- phy@edge[!is.na(phy@edge[, 1]), ]
    ## Sort edges -- ensures that starting order of edge matrix doesn't
    ## affect the order of reordered trees
    edge <- edge[order(edge[, 2]), ]

    # recast order argument as integer to pass to C
    if(order == 'postorder') {
        iOrder <- 0L
    } else if(order == 'preorder') {
        iOrder <- 1L
    } else {stop(paste("Method for", order, "not implemented"))}

    if (!hasPoly(phy) & !hasSingle(phy)) {
        # method 1: faster, but only works if all internal nodes have
        # exactly two children (true binary tree)

        # extract nodes, separating descendants into left (first
        # encountered) and right (second encountered) for each ancestor
        isFirst <- !duplicated(edge[, 1])
        ancestor <- as.integer(edge[isFirst, 1])
        left <- as.integer(edge[isFirst, 2])
        right <- as.integer(edge[!isFirst, 2])[match(ancestor,
            edge[!isFirst, 1])]
        descendantNew <- rep(0L, nEdges(phy))
        root <- as.integer(rootNode(phy))
        nEdge <- as.integer(length(ancestor))

        descendantReord <- .C("reorderBinary", descendantNew, root,
            ancestor, left, right, nEdge, iOrder)[[1]]

    } else {
        # method 2: not as fast, but robust to singletons and polytomies

        # extract ancestors and descendants
        ancestor <- as.integer(edge[,1])
        descendant <- as.integer(edge[,2])
        descendantNew <- rep(0L, nEdges(phy))
        root <- as.integer(rootNode(phy))
        nEdge <- as.integer(nrow(edge))

        descendantReord <- .C("reorderRobust", descendantNew, root,
            ancestor, descendant, nEdge, iOrder)[[1]]

    }

    ## Original pure R implementation of the above:
    #### recursive functions are placed first and calls to those functions below
    ##postOrder <- function(node) {
    ##    ## this function returns a vector of nodes in the post order traversal
    ##    ## get the descendants
    ##    traversal <- NULL
    ##    ## edge -- defined above, outside this function
    ##    ## extensive testing found this loop to be faster than apply() etc.
    ##    for(i in edge[edge[, 1] == node, 2]) {
    ##        traversal <- c(traversal, postOrder(i))
    ##    }
    ##    c(traversal, node)
    ##}
    ##preOrder  <- function(node) {
    ##    ## see expanded code in comments of postOrder()
    ##    ## only difference here is that we record current node, then descendants
    ##    traversal <- NULL
    ##    for(i in edge[edge[, 1] == node, 2]) {
    ##        traversal <- c(traversal, preOrder(i))
    ##    }
    ##    c(node, traversal)
    ##}
    ##if(order == 'postorder') {
    ##    descendantReord <- postOrder(rootNode(phy))
    ##} else if(order == 'preorder') {
    ##    descendantReord <- preOrder(rootNode(phy))
    ##} else {stop(paste("Method for", order, "not implemented"))}

    ## match the new node order to the old order to get an index
    index <- match(descendantReord, phy@edge[, 2])

}

setMethod("reorder", signature(x = 'phylo4'),
    function(x, order = c('preorder', 'postorder')) {
    ## call orderIndex and use that index to order edges, labels and lengths
    order   <- match.arg(order)
    index   <- orderIndex(x, order)
    x@order <- order
    x@edge  <- x@edge[index, ]
    if(hasEdgeLabels(x)) { x@edge.label  <- x@edge.label[index] }
    if(hasEdgeLength(x)) { x@edge.length <- x@edge.length[index] }
    x
})


