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

setMethod("nTips", signature(x="phylo4"), function(x) {
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
setMethod("nTips", signature(x="phylo"),
 function(x) {
     Ntip(x)
})

#########################################################
### Node accessors
#########################################################

setMethod("nNodes", signature(x="phylo4"), function(x) {
    E <- edges(x, drop.root=TRUE)
    if(nrow(E) == 0) {
        return(0)
    } else {
        return(length(unique(E[, 1])))
    }
})

setMethod("nodeType", signature(x="phylo4"),
 function(x) {
    if(nTips(x) == 0)
        return(NULL)
    else {
        ## strip out the root ancestor
        nodesVect <- as.vector(edges(x))
        nodesVect <- nodesVect[nodesVect != 0]
        ## get a sorted list of the unique nodes 
        listNodes <- sort(unique(nodesVect))
        t <- rep("internal", length(listNodes)) # FM: internal is default (I think it's safer)
        names(t) <- listNodes

        ## node number of real internal nodes
        iN <- names(table(edges(x)[,1]))
        ## node number that are not internal nodes (ie that are tips)
        tN <- names(t)[!names(t) %in% iN]
        t[tN] <- "tip"

        ## if the tree is rooted
        if(isRooted(x)) t[rootNode(x)] <- "root"

        return(t)
    }
})

# return node IDs (or a subset thereof) in ascending order
setMethod("nodeId", signature(x="phylo4"),
 function(x, type=c("all",
    "tip","internal","root")) {

     type <- match.arg(type)
     E <- edges(x)

     ## Note: this implementation will still work even if tips are not
     ## 1:nTips and nodes are not (nTips+1):nNodes
     nid <- switch(type,
         ## all nodes appear at least once in the edge matrix
         all = unique(as.vector(E)[as.vector(E) != 0]),
         ## tips are nodes that do not appear in the ancestor column
         tip = setdiff(E[, 2], E[, 1]),
         ## internals are nodes that *do* appear in the ancestor column
         internal = unique(E[E[, 1] != 0, 1]),
         ## roots are nodes that have NA as ancestor
         root = if (!isRooted(x)) NA else unname(E[E[, 1] == 0, 2]))

     return(sort(nid))

})



#########################################################
### Edge accessors
#########################################################

setMethod("nEdges", signature(x="phylo4"),
 function(x) {
    nrow(x@edge)
})

# return edge matrix in its current order
setMethod("edges", signature(x="phylo4"),
 function(x, drop.root=FALSE, ...) {
     e <- x@edge
     if (drop.root) e <- e[e[, 1] != 0, ]
     e
})

setMethod("edgeOrder", signature(x="phylo4"),
 function(x, ...) {
    x@order
})

# return edge IDs (or a subset thereof) in edge matrix order
setMethod("edgeId", signature(x="phylo4"),
 function(x, type=c("all", "tip",
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
        isRoot <- edge[, 1] == 0
        edge <- edge[isRoot, , drop=FALSE]
    } # else just use complete edge matrix if type is "all"
    id <- paste(edge[, 1], edge[, 2], sep="-")
    return(id)
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
        x@edge.length <- numeric()
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

#########################################################
### Root accessors
#########################################################

setMethod("isRooted", signature(x="phylo4"),
 function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x) == 0) return(FALSE)
    any(edges(x)[, 1] == 0)
})

setMethod("rootNode", signature(x="phylo4"),
 function(x) {
    if (!isRooted(x))
        return(NA)
    unname(edges(x)[which(edges(x)[, 1] == 0), 2])
})

setReplaceMethod("rootNode", signature(x="phylo4"),
 function(x, value) {
    stop("Root node replacement not implemented yet")
})

#########################################################
### Label accessors
#########################################################

## return labels in increasing node order
setMethod("labels", signature(object="phylo4"),
  function(object, type = c("all", "tip", "internal")) {
    type <- match.arg(type)
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    label <- object@label
    id <- nodeId(object, type)
    lbl <- label[match(id, names(label))]
    # reassign names b/c any unmatched will be NA (could instead assign
    # names only to the unmatched ones, but this seems simpler)
    names(lbl) <- id
    return(lbl)
})

setReplaceMethod("labels",
                 signature(x="phylo4", type="ANY",
                           use.names="ANY", value="character"),
   function(x, type = c("all", "tip", "internal"),
            use.names, ..., value) {

       ## Default options
       if(missing(type))
           type <- "all"
       if (missing(use.names))
           use.names <- FALSE

       type <- match.arg(type)

       ## generate new labels of the desired type
       new.label <- .createLabels(value, nTips(x), nNodes(x), use.names,
           type=type)

       ## replace existing labels and add new ones as needed
       old.label <- x@label
       old.index <- match(names(new.label), names(old.label))
       isNew <- is.na(old.index)
       old.label[old.index[!isNew]] <- new.label[!isNew]
       updated.label <- c(old.label, new.label[isNew])

       ## for efficiency, drop any NA labels
       x@label <- updated.label[!is.na(updated.label)]

       if(is.character(checkval <- checkPhylo4(x)))
           stop(checkval)
       else
           return(x)
   })


### Node Labels
setMethod("hasNodeLabels", signature(x="phylo4"),
 function(x) {
    !all(is.na(nodeLabels(x)))
})

setMethod("nodeLabels", signature(x="phylo4"),
 function(x) {
    labels(x, type="internal")
})

setReplaceMethod("nodeLabels", signature(x="phylo4", value="character"),
  function(x, ..., value) {
      labels(x, type="internal", ...) <- value
      if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
      x
  })

### Tip labels
setMethod("tipLabels", signature(x="phylo4"),
 function(x) {
    labels(x, type="tip")
    })

setReplaceMethod("tipLabels", signature(x="phylo4", value="character"),
  function(x, ...,  value) {
      labels(x, type="tip", ...) <- value
      if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
      return(x)
  })


### Edge labels
setMethod("hasEdgeLabels", signature(x="phylo4"),
 function(x) {
    !all(is.na(x@edge.label))
})

# return edge labels in order by edgeIds (same order as edge matrix)
setMethod("edgeLabels", signature(x="phylo4"),
  function(x) {
    ## [JR: below, using match for ordering rather than direct character
    ## indexing b/c the latter is slow for vectors of a certain size]
    id <- edgeId(x, "all")
    lbl <- x@edge.label[match(id, names(x@edge.label))]
    names(lbl) <- id
    return(lbl)
})

setReplaceMethod("edgeLabels", signature(x="phylo4", value="character"),
  function(x, ..., value) {
    lbl <- .createEdge(value, x@edge, type="labels")
    x@edge.label <- lbl[!is.na(lbl)]
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    x
  })


#########################################################
### Displaying phylo4: print, show, head, tail and summary
#########################################################

### print
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

### summary
setMethod("summary", signature(object="phylo4"),
  function(object, quiet=FALSE) {

      res <- list()

      ## build the result object
      res$name <- deparse(substitute(object, sys.frame(-1)))
      res$nb.tips <- nTips(object)
      res$nb.nodes <- nNodes(object)

      if(hasEdgeLength(object)) {
          edge.length <- edgeLength(object)
          res$mean.el <- mean(edge.length, na.rm=TRUE)
          res$var.el <- var(edge.length, na.rm=TRUE)
          if (isRooted(object) && is.na(edgeLength(object, rootNode(object)))) {
              root.index <- match(edgeId(object, "root"), names(edge.length))
              res$sumry.el <- summary(edge.length[-root.index])
          } else {
              res$sumry.el <- summary(edge.length)
          }
      }

      ## check for polytomies
      if (hasPoly(object)) {
          E <- edges(object)
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
          names(res$degree) <- nodeLabels(object)
          res$polytomy <- res$polytomy[idx]
          names(res$polytomy) <- nodeLabels(object)
      }

      ## model info
      res$loglik <- attr(object, "loglik")
      res$para <- attr(object, "para")
      res$xi <- attr(object, "xi")

      ## if quiet, stop here
      if(quiet) return(invisible(res))

      ## now, print to screen is !quiet
      cat("\n Phylogenetic tree :", res$name, "\n\n")
      cat(" Number of tips    :", res$nb.tips, "\n")
      cat(" Number of nodes   :", res$nb.nodes, "\n")
      ## cat("  ")
      if(hasEdgeLength(object)) {
          cat(" Branch lengths:\n")
          cat("        mean         :", res$mean.el, "\n")
          cat("        variance     :", res$var.el, "\n")
          cat("        distribution :\n")
          print(res$sumry.el)
      }
      else {
          cat(" Branch lengths    : No branch lengths.\n")
      }
      if (hasPoly(object)) {
          cat("\nDegree of the nodes  :\n")
          print(res$degree)
          cat("\n")
          cat("Types of polytomy:\n")
          print(res$polytomy)
          cat("\n")
      }

      if (!is.null(attr(object, "loglik"))) {
          cat("Phylogeny estimated by maximum likelihood.\n")
          cat("  log-likelihood:", attr(object, "loglik"), "\n\n")
          npart <- length(attr(object, "para"))
          for (i in 1:npart) {
              cat("partition ", i, ":\n", sep = "")
              print(attr(object, "para")[[i]])
              if (i == 1)
                  next
              else cat("  contrast parameter (xi):", attr(object,"xi")[i - 1], "\n")
        }
      }
      return(invisible(res))

  }) # end setMethod summary phylo4




#########################################################
### Ordering
#########################################################

orderIndex <- function(x, order=c("preorder", "postorder")) {

    browser()
    order <- match.arg(order)
    if(!isRooted(x)){
        stop("Tree must be rooted to reorder")
    }
    ## get a root node free edge matrix
    edge <- edges(x, drop.root=TRUE)
    ## Sort edges -- ensures that starting order of edge matrix doesn't
    ## affect the order of reordered trees
    edge <- edge[order(edge[, 2]), ]

    # recast order argument as integer to pass to C
    if(order == 'postorder') {
        iOrder <- 0L
    } else if(order == 'preorder') {
        iOrder <- 1L
    } else {stop(paste("Method for", order, "not implemented"))}

    if (!hasPoly(x) & !hasSingle(x)) {
        # method 1: faster, but only works if all internal nodes have
        # exactly two children (true binary tree)

        # extract nodes, separating descendants into left (first
        # encountered) and right (second encountered) for each ancestor
        isFirst <- !duplicated(edge[, 1])
        ancestor <- as.integer(edge[isFirst, 1])
        left <- as.integer(edge[isFirst, 2])
        right <- as.integer(edge[!isFirst, 2])[match(ancestor,
            edge[!isFirst, 1])]
        descendantNew <- rep(0L, nEdges(x))
        root <- as.integer(rootNode(x))
        nEdge <- as.integer(length(ancestor))

        descendantReord <- .C("reorderBinary", descendantNew, root,
            ancestor, left, right, nEdge, iOrder)[[1]]

    } else {
        # method 2: not as fast, but robust to singletons and polytomies

        # extract ancestors and descendants
        ancestor <- as.integer(edge[,1])
        descendant <- as.integer(edge[,2])
        descendantNew <- rep(0L, nEdges(x))
        root <- as.integer(rootNode(x))
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
    ##    descendantReord <- postOrder(rootNode(x))
    ##} else if(order == 'preorder') {
    ##    descendantReord <- preOrder(rootNode(x))
    ##} else {stop(paste("Method for", order, "not implemented"))}

    ## match the new node order to the old order to get an index
    index <- match(descendantReord, edges(x)[, 2])

}

setMethod("reorder", signature(x="phylo4"),
 function(x, order=c("preorder", "postorder")) {
    ## call orderIndex and use that index to order edges, labels and lengths
    order   <- match.arg(order)
    index   <- orderIndex(x, order)
    x@order <- order
    x@edge  <- edges(x)[index, ]
    if(hasEdgeLabels(x)) {
        x@edge.label  <- x@edge.label[index]
    }
    if(hasEdgeLength(x)) {
        x@edge.length <- x@edge.length[index]
    }
    x
})


