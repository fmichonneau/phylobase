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
###  3.4. hasEdgeLength()
###  3.5. edgeLength()
###  3.6. edgeLength() <-
###  3.7. sumEdgeLength()

### 4. Root accessors
###  4.1. isRooted()
###  4.2. rootNode()

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

nTips <- function(x,...)  { }  ## mask ape::nTips
setGeneric("nTips", function(x,...) {
    standardGeneric("nTips")
})

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

setMethod("nodeId", "phylo4", function(x, type=c("internal","tip","allnode")) {
  type <- match.arg(type)
  tipNid <- x@edge[x@edge[,2]<=nTips(x),2]
  allNid <- unique(as.vector(x@edge))
  intNid <- allNid[! allNid %in% tipNid]
  nid <- switch(type,
                internal = intNid,
                tip = tipNid,
                allnode = allNid)
  return(nid[!is.na(nid)])
})



#########################################################
### Edge accessors
#########################################################

setMethod("nEdges", "phylo4", function(x) {
    nrow(x@edge)
})

setMethod("edges", "phylo4", function(x, order, drop.root=FALSE, ...) {
  e <- x@edge
  if (drop.root) e <- e[!is.na(e[,1]),]
  e
})

setMethod("edgeOrder", "phylo4", function(x, ...) {
    x@order
})

setMethod("hasEdgeLength","phylo4", function(x) {
    !all(is.na(x@edge.length))
})

setMethod("edgeLength", "phylo4", function(x, node) {
    if (!hasEdgeLength(x))
        NULL
    else {
      if (missing(node))
          return(x@edge.length)
      else {
          n <- getNode(x, node)
          return(x@edge.length[match(n, x@edge[,2])])
      }
    }
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

setMethod("labels", "phylo4", function(object, type = c("tip",
    "internal", "allnode"), ...) {
    type <- match.arg(type)
    switch(type,
            tip = object@tip.label[as.character(nodeId(object, "tip"))],
            internal = {
                if (hasNodeLabels(object)) {
                    object@node.label
                }
                else
                {
                    ## FIXME? should this return object@node.label
                    return(character(0))
                }
            },
            allnode = {
                c(object@tip.label, object@node.label)
              }
            )
})

setReplaceMethod("labels",
                 signature(object="phylo4", value="character"),
   function(object, type = c("tip", "internal", "allnode"),
            use.names=FALSE, ..., value) {

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
      object
  })

### Tip labels
setMethod("tipLabels", "phylo4", function(object) {
    labels(object, type="tip")
    })

setReplaceMethod("tipLabels", signature(object="phylo4", value="character"),
  function(object, ...,  value) {
      labels(object, type="tip", ...) <- value
      return(object)
  })


### Edge labels
setMethod("hasEdgeLabels", "phylo4", function(x) {
    length(x@edge.label) > 0
})

setMethod("edgeLabels", signature(x = "phylo4"), function(x) {
    x@edge.label
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
printphylo4 <- function(x, printall=TRUE) {
    if(!nrow(x@edge)) {
        msg <- paste("Empty \'", class(x), "\' object\n", sep="")
        cat(msg)
    }
    else {
        if (printall)
            print(as(x, 'data.frame'))
        else print(head(as(x, 'data.frame')))
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
        res$sumry.el <- summary(x@edge.length)[-4]
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
    ## get an root node free edge matrix
    ## R scoping allows us to call this variable in
    ## the postOrder() func defined above
    edge <- phy@edge[!is.na(phy@edge[, 1]), ]
    ## Sort edges -- ensures that starting order of edge matrix doesn't
    ## affect the order of reordered trees
    edge <- edge[order(edge[, 2]), ]

    ## recursive functions are placed first and calls to those functions below
    postOrder <- function(node) {
        ## this function returns a vector of nodes in the post order traversal
        ## get the descendants
        traversal <- NULL
        ## edge -- defined above, outside this function
        ## extensive testing found this loop to be faster than apply() etc.
        for(i in edge[edge[, 1] == node, 2]) {
            traversal <- c(traversal, postOrder(i))
        }
        c(traversal, node)
    }
    preOrder  <- function(node) {
        ## see expanded code in comments of postOrder()
        ## only difference here is that we record current node, then descendants
        traversal <- NULL
        for(i in edge[edge[, 1] == node, 2]) {
            traversal <- c(traversal, preOrder(i))
        }
        c(node, traversal)
    }

    if(order == 'postorder') {
        ## match the new node order to the old order to get an index
        index <- match(postOrder(rootNode(phy)), phy@edge[, 2])

    } else if(order == 'preorder') {
        ## match the new node order to the old order to get an index
        index <- match(preOrder(rootNode(phy)), phy@edge[, 2])

    } else {stop(paste("Method for", order, "not implemented"))}
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


