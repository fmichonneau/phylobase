#######################################################
## Importing from ape
setAs("phylo", "phylo4", function(from, to) {
    ## fixme SWK kludgy fix may not work well with unrooted trees
    ## TODO should we also attempt to get order information?
    if (is.rooted(from)) {
        tip.idx <- 1:nTips(from)
        int.idx <- (nTips(from)+1):dim(from$edge)[1]
        root.node <- as.numeric(setdiff(unique(from$edge[,1]), unique(from$edge[,2])))

        from$edge <- rbind(from$edge[tip.idx,],c(NA,root.node),from$edge[int.idx,])
        if (!is.null(from$edge.length)) {
            if (is.null(from$root.edge)) {
                from$edge.length <- c(from$edge.length[tip.idx],as.numeric(NA),from$edge.length[int.idx])
            }
            else {
                from$edge.length <- c(from$edge.length[tip.idx],from$root.edge,from$edge.length[int.idx])
            }
        }
        if (!is.null(from$edge.label)) {
            from$edge.label <- c(from$edge.label[tip.idx],NA,from$edge.label[int.idx])
        }
    }
    newobj <- phylo4(from$edge, from$edge.length, from$tip.label,
                     node.label = from$node.label, edge.label = from$edge.label)
    attribs <- attributes(from)
    attribs$names <- NULL
    knownattr <- c("logLik", "order", "origin", "para", "xi")
    known <- names(attribs)[names(attribs) %in% knownattr]
    unknown <- names(attribs)[!names(attribs) %in% c(knownattr, "class", "names")]
    if (length(unknown) > 0) {
        warning(paste("unknown attributes ignored: ", unknown, collapse = " "))
    }
    for (i in known) attr(newobj, i) <- attr(from, i)
    newobj
})

setAs("phylo", "phylo4d", function(from, to) {
    phylo4d(as(from, "phylo4"), tip.data = data.frame())
})


setAs("multiPhylo", "multiPhylo4", function(from, to) {
    trNm <- names(from)
    if(is.null(trNm)) trNm <- character(0)
    newobj <- new("multiPhylo4", phylolist = lapply(from, function(x)
                                 as(x, "phylo4")),
                  tree.names = trNm)
    newobj
})

#######################################################
## Exporting to ape


## BMB: adding an explicit as method, and the warning,
##  here is a very bad idea, because
##   even implicit conversions from phylo4d to phylo4 (e.g.
##  to use inhertied methods) will produce the warning

## setAs("phylo4d", "phylo4",function(from,to) {
##   warning("losing data while coercing phylo4d to phylo")
##   phylo4(from@edge, from@edge.length, from@tip.label,
##         from@node.label,from@edge.label,from@order)
## })

setAs("phylo4", "phylo", function(from, to) {
    if (inherits(from, "phylo4d"))
        warning("losing data while coercing phylo4d to phylo")
    brlen <- from@edge.length
    ## rootnode is only node with no ancestor
    rootpos <- which(is.na(from@edge[, 1]))
    if (isRooted(from)) brlen <- brlen[-rootpos]
    edgemat <- unname(from@edge[-rootpos, ])
    y <- list(edge = edgemat,
            Nnode = from@Nnode,
            tip.label = from@tip.label,
            edge.length = brlen,
            node.label = from@node.label)
    class(y) <- "phylo"
    if (from@order != 'unknown') {
        ## TODO postorder != pruningwise -- though quite similar
        attr(y, 'order') <- switch(from@order, postorder = 'unknown',
                                      preorder  = 'cladewise')
    }
    if (length(y$edge.length) == 0)
        y$edge.length <- NULL
    if (length(y$node.label) == 0)
        y$node.label <- NULL
    if (isRooted(from)) {
        root.edge <- brlen[nodeId(from, "all") == rootNode(from)]
        if (!is.na(root.edge)) y$root.edge <- root.edge
    }
    y
})

## BMB: redundant????
## setAs("phylo4d", "phylo", function(from, to) {
##     y <- list(edge = from@edge, edge.length = from@edge.length,
##         Nnode = from@Nnode, tip.label = from@tip.label)
##     class(y) <- "phylo"
##     if (length(y$edge.length) == 0)
##         y$edge.length <- NULL
##     if (length(y$node.label) == 0)
##         y$node.label <- NULL
##     #if (!is.na(from@root.edge))
##     #    y$root.edge <- from@root.edge
##    warning("losing data while coercing phylo4d to phylo")
##    y
##})

setAs("multiPhylo4", "multiPhylo", function(from, to) {
    y <- lapply(from@phylolist, function(x) as(x, "phylo"))
    names(y) <- from@tree.names
    if (nrow(from@tip.data) > 0)
        warning("discarded tip data")
    class(y) <- "multiPhylo"
    y
})

#######################################################
## Exporting to ade4
setAs("phylo4", "phylog", function(from, to) {
    if (!require(ade4))
        stop("the ade4 package is required")
    x <- as(from, "phylo")
    xstring <- write.tree(x, file = "")
    newick2phylog(xstring)
})

#######################################################
## Exporting to dataframe
setAs(from = "phylo4", to = "data.frame", def = function(from) {

    ## Check the phylo4
    if (is.character(checkval <- checkPhylo4(from)))
        stop(checkval)
    x <- from

    node <- nodeId(x, "all")
    ancestr <- ancestor(x, node)
    ndType <- nodeType(x)
    intNode <- names(ndType[ndType == "internal"])
    tip <- names(ndType[ndType == "tip"])

    E <- data.frame(node, ancestr)
    ## !! in phylobase, the order is node-ancestors whereas in ape it's
    ## ancestor-node
    nmE <- paste(E[,2], E[,1], sep="-")

    if (hasEdgeLength(x)) {
        edge.length <- edgeLength(x)[match(nmE, names(x@edge.length))]
    }
    else {
        edge.length <- rep(NA, nEdges(x))
    }

    label <- labels(x,which="all")[node]

    d <- data.frame(label, node, ancestr, edge.length, node.type=ndType[node],
                    row.names=node)
    d$label <- as.character(d$label)
    d
})

setAs(from = "phylo4d", to = "data.frame", function(from) {

    tree <- extractTree(from)
    ## Convert to data.frame
    tDf <- as(tree, "data.frame")

    dat <- tdata(from, "allnode", label.type="column") # get data

    ## reorder data to edge matrix order, drop labels (first column)
    if(nrow(dat) > 0 && ncol(dat) > 1) {
         dat <- dat[match(rownames(tDf), rownames(dat)), ]
         tdat <- cbind(tDf, dat[ ,-1 , drop=FALSE])
     }
     else {
        tdat <- tDf
        cat("No data associated with the tree\n")
     }
    tdat
})
