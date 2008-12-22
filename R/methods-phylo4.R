## accessor functions for all internal bits
## HORRIBLE KLUGE
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

setMethod("nNodes", "phylo4", function(x) {
    x@Nnode
})

setMethod("nEdges", "phylo4", function(x) {
    nrow(x@edge)
})

setMethod("edges", "phylo4", function(x, order, ...) {
    x@edge
})

setMethod("edgeOrder", "phylo4", function(x, ...) {
    x@order
})


setMethod("isRooted","phylo4", function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x) == 0) return(FALSE)
    any(is.na(edges(x)[,1]))
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


setMethod("rootNode", "phylo4", function(x) {
    if (!isRooted(x))
        return(NA)
    unname(edges(x)[which(is.na(edges(x)[,1])),2])
})

setReplaceMethod("rootNode", "phylo4", function(x, value) {
    stop("Root node replacement not implemented yet")
})

setMethod("edgeLength", "phylo4", function(x,which) {
    if (!hasEdgeLength(x))
        NULL
    else {
      if (missing(which)) return(x@edge.length)
      n <- getNode(x,which)
      return(x@edge.length[n])
    }
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

setMethod("hasNodeLabels", "phylo4", function(x) {
    length(x@node.label) > 0
})

setMethod("hasEdgeLabels", "phylo4", function(x) {
    length(x@edge.label) > 0
})

setMethod("labels", "phylo4", function(object, which = c("tip",
    "node", "allnode"), ...) {
    which <- match.arg(which)
    switch(which,
            tip = object@tip.label[order(nodeId(object,"tip"))],
            node = {
                if (hasNodeLabels(object)) {
                    object@node.label[order(nodeId(object))]
                }
                else
                {
                    return(character(0))
                }
            }
            ,
            allnode = {
                if (hasNodeLabels(object)) {
                    nl <- object@node.label
                }
                else
                {
                    nl <- rep(NA,nNodes(object))
                }
                ## lorder <- match(object@edge[,2],
                ## c(nodeId(object,"tip"),nodeId(object)))
                lorder <- order(c(nodeId(object,"tip"),nodeId(object)))
                c(object@tip.label,nl)[lorder]
              }
            )
})

setMethod("nodeLabels", "phylo4", function(x) {
    x@node.label
})

setMethod("nodeId", "phylo4", function(x,which=c("internal","tip","all")) {
  which <- match.arg(which)
  switch(which,
         internal=x@edge[x@edge[,2]>nTips(x),2],
         tip = x@edge[x@edge[,2]<=nTips(x),2],
         all = x@edge[,2])
})

setReplaceMethod("nodeLabels", "phylo4",
                 function(object, ..., value) {
                   ## FIXME: test length!
                   if (length(value)!=nNodes(object))
                     stop("label vector must have as many elements as number of internal nodes")
                   object@node.label <- value
                   object
                 })

setMethod("edgeLabels", "phylo4", function(x) {
    x@edge.label
})

setReplaceMethod("edgeLabels", "phylo4", function(object, ...,
    value) {
    object@edge.label <- value
    object
})

## hack to allow access with $
setMethod("$", "phylo4", function(x, name) {
    switch(name, edge.length = if (!hasEdgeLength(x))
        NULL
    else x@edge.length, node.label = if (!hasNodeLabels(x))
        NULL
    else x@node.label, attr(x, name))
})

## FIXME: implement more checks on this!!
setReplaceMethod("$", "phylo4", function(x, name, value) {
    slot(x, name, check = TRUE) <- value
    return(x)
})


printphylo4 <- function(x, printall = TRUE){
    if (printall)
      print(as(x, 'data.frame'))
    else print(head(as(x, 'data.frame')))
}
## hack for print/show
## from http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1363.html
#setMethod("print", "phylo4", printphylo)
#setMethod("show", "phylo4", function(object) printphylo(object))
setMethod("print", "phylo4", printphylo4)
setMethod("show", "phylo4", function(object) printphylo4(object))
##
# Alternative print method for phylo4, showing the contents of the tree data.
##  Not sure if it works for unrooted trees

printphylo <- function (x,printlen=6,...) {
    printlen <- max(1,printlen)
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    nb.edge <- length(x$edge.label)
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and",
              nb.node, "internal nodes\n"))

    ## print tip labels
    cat("\nTip labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x$tip.label[1:printlen], collapse = ", "),
                  ", ...\n", sep = ""))
    } else print(x$tip.label)

    ## print node labels
    cat("\nNode labels:\n")
    if (nb.node > printlen) {
        cat(paste("\t", paste(x$node.label[1:printlen], collapse = ", "),
                  ", ...\n", sep = ""))
    } else print(x$node.label)

    ## print edge labels
    cat("\nEdge labels:\n")
    if (nb.edge > printlen) {
        cat(paste("\t", paste(x$edge.label[1:printlen], collapse = ", "),
                  ", ...\n", sep = ""))
    } else print(x$edge.label)

    ## slots
    ##     cat("\nSlots:\n")
    ##     cat(paste("@", names(x)[1:4], sep=""),sep="\t")
    ##     cat("\n")
    ##     cat(paste("@", names(x)[5:7], sep=""),sep="\t")
    ##     cat("\n")

    rlab <- if (isRooted(x)) "Rooted"  else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (hasEdgeLength(x))
      "includes branch lengths"
    else       "no branch lengths"
    cat(blen, "\n\n", sep = "")
}

#################
## summary phylo4
#################
setMethod("summary","phylo4", function (object, quiet=FALSE) {
    x <- object
    res <- list()

    ## build the result object
    res$name <- deparse(substitute(object, sys.frame(-1)))
    res$nb.tips <- length(x$tip.label)
    res$nb.nodes <- x$Nnode

    if(!is.null(x$edge.length)){
        res$mean.el <- mean(x$edge.length, na.rm=TRUE)
        res$var.el <- var(x$edge.length, na.rm=TRUE)
        res$sumry.el <- summary(x$edge.length)[-4]
    } else {
        res$mean.el <- NULL
        res$var.el <- NULL
        res$sumry.el <- NULL
    }

    ## check for polytomies
    if (any(tabulate(na.omit(edges(object)[,1]))>2)){ # if there are polytomies
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
    if(is.null(x$edge.length)) {
        cat(" Branch lengths    : No branch lengths.\n")
    } else {
        cat(" Branch lengths:\n")
        cat("        mean         :", res$mean.el, "\n")
        cat("        variance     :", res$var.el, "\n")
        cat("        distribution :\n")
        print(res$sumry.el)
    }
    if(hasPoly(x)){
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

setMethod("names", signature(x = "phylo4"), function(x){
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})

setMethod("hasEdgeLength","phylo4", function(x) {
    length(x@edge.length)>0
})

setReplaceMethod("labels","phylo4", function(object,...,value) {
    if (length(value) != length(object@tip.label))
        stop("Number of tip labels does not match number of tips.")
    object@tip.label <- value
    object
})

orderIndex <- function(phy, order = c('preorder', 'postorder')) {
    ## recursive functions are placed first and calls to those functions below
    postOrder <- function(node) {
        ## this function returns a list of nodes in the post order traversal
        ## get the descendants 
        ## dec <- edge[, 1] == node
        ## print(dec)
        ## recursive call to get the descendants of the descendants
        ## out <- mapply(postie, edge[dec, 2])
        ## return the descendants with the node after
        ## return(c(unlist(out), node))
        ## slight performance benefit to the one liner
        return(c(unlist(mapply(postOrder, edge[edge[, 1] == node, 2])), node))
    }
    preOrder  <- function(node) {
        ## see expanded code in comments of postOrder
        ## only difference here is that we record current node, then descendants
        return(c(node, unlist(mapply(preOrder, edge[edge[, 1] == node, 2]))))
    }
    
    if(order == 'postorder') {
        ## get an root node free edge matrix 
        edge <- phy@edge[!is.na(phy@edge[, 1]), ]
        ## match the new node order to the old order to get an index
        index <- match(postOrder(rootNode(phy)), phy@edge[, 2])
    
    } else if(order == 'preorder') {
        ## get an root node free edge matrix 
        edge <- phy@edge[!is.na(phy@edge[, 1]), ]
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

setMethod("head",signature(x = 'phylo4'),
          function(x,n=20) {
            head(as(x,"data.frame"),n=n)
          })

setMethod("tail",signature(x = 'phylo4'),
          function(x,n=20) {
            tail(as(x,"data.frame"),n=n)
          })
