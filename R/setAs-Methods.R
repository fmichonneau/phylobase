#######################################################
## Importing from ape
setAs("phylo", "phylo4", function(from, to) {
    ## fixme SWK kludgy fix may not work well with unrooted trees
    ## TODO should we also attempt to get order information?
    if (is.rooted(from)) {
        tip.idx <- 1:nTips(from)
        if (nTips(from) < nrow(from$edge)) {
            int.idx <- (nTips(from)+1):dim(from$edge)[1]
        } else {
            int.idx <- NULL
        }
        root.node <- as.numeric(setdiff(unique(from$edge[,1]), unique(from$edge[,2])))

        from$edge <- rbind(from$edge[tip.idx,],c(0,root.node),from$edge[int.idx,])
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
    oldorder <- attr(from,"order")
    neworder <- if (is.null(oldorder)) { "unknown" } else
    if (!oldorder %in% phylo4_orderings) {
      stop("unknown ordering '",oldorder,"' in ape object")
    } else if (oldorder=="cladewise") "preorder"
    else oldorder
    attr(from,"order") <- NULL
    newobj <- phylo4(from$edge, from$edge.length, from$tip.label,
                     node.label = from$node.label,
                     edge.label = from$edge.label,
                     order = neworder)
    attribs <- attributes(from)
    attribs$names <- NULL
    knownattr <- c("logLik", "origin", "para", "xi")
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
##  to use inherited methods) will produce the warning

## setAs("phylo4d", "phylo4",function(from,to) {
##   warning("losing data while coercing phylo4d to phylo")
##   phylo4(from@edge, from@edge.length, from@tip.label,
##         from@node.label,from@edge.label,from@order)
## })

setAs("phylo4", "phylo", function(from, to) {

    if(is.character(checkval <- checkPhylo4(from))) {
      stop(checkval)
    }

    if (inherits(from, "phylo4d"))
        warning("losing data while coercing phylo4d to phylo")

    phy <- list()

    ## Edge matrix (dropping root edge if it exists)
    edgemat <- unname(edges(from, drop.root=TRUE))
    storage.mode(edgemat) <- "integer"
    phy$edge <- edgemat

    ## Edge lengths
    if(hasEdgeLength(from)) {
        edge.length <- edgeLength(from)
        if(isRooted(from)) {
            iRoot <- match(edgeId(from, "root"), names(edge.length))
            phy$edge.length <- unname(edge.length[-iRoot])
        }
        else {
            phy$edge.length <- unname(edge.length)
        }
    }

    ## Tip labels
    phy$tip.label <- unname(tipLabels(from))

    ## nNodes
    phy$Nnode <- as.integer(nNodes(from))

    ## Node labels
    if(hasNodeLabels(from)) {
        phy$node.label <- unname(nodeLabels(from))
    }

    ## Root edge
    if(isRooted(from) && hasEdgeLength(from)) {
        root.edge <- unname(edgeLength(from,rootNode(from)))
        if(!is.na(root.edge)) {
            phy$root.edge <- root.edge
        }
    }

    ## Converting to class phylo
    class(phy) <- "phylo"

    ## Tree order
    ## TODO postorder != pruningwise -- though quite similar
    if (edgeOrder(from) == "unknown") {
        warning("trees with unknown order may be",
                " unsafe in ape")
    }
    else {
        attr(phy, "order") <- switch(edgeOrder(from),
                                     postorder = "unknown",
                                     preorder = "cladewise",
                                     pruningwise = "pruningwise")
    }
    phy
})


## BMB: redundant????
## setAs("phylo4d", "phylo", function(from, to) {
##     y <- list(edge = from@edge, edge.length = from@edge.length,
##         Nnode = nNodes(from), tip.label = from@tip.label)
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

.phylo4ToDataFrame <- function(from, edgeOrder=c("pretty", "real")) {

    edgeOrder <- match.arg(edgeOrder)

    ## Check the phylo4
    if (is.character(checkval <- checkPhylo4(from)))
        stop(checkval)

    ## The order of 'node' defines the order of all other elements
    if (edgeOrder == "pretty") {
        node <- nodeId(from, "all")
        ancestr <- ancestor(from, node)

        # ancestor returns an NA, replace this w/ 0 to construct names correctly
        ancestr[is.na(ancestr)] <- as.integer(0)
    } else {
        E <- edges(from)
        node <- E[, 2]
        ancestr <- E[, 1]
    }

    ## extract and reorder (as needed) other object slots
    nmE <- paste(ancestr, node, sep="-")
    edge.length <- edgeLength(from)
    edge.length <- edge.length[match(nmE, names(edge.length))]

    ndType <- nodeType(from)
    ndType <- ndType[match(node, names(ndType))]
    label <- labels(from, type="all")
    label <- label[match(node, names(label))]

    tDf <- data.frame(label, node, ancestor=ancestr, edge.length,
                    node.type=ndType, row.names=node)
    tDf$label <- as.character(tDf$label)

    if (class(from) == "phylo4d") {
        dat <- tdata(from, "allnode", label.type="column") # get data

        ## reorder data to edge matrix order, drop labels (first column)
        if(nrow(dat) > 0 && ncol(dat) > 1) {
            dat <- dat[match(rownames(tDf), rownames(dat)), ]
            tDf <- cbind(tDf, dat[ ,-1 , drop=FALSE])
        }
        else {
            cat("No data associated with the tree\n")
       }
    }
    tDf
}

setAs(from = "phylo4", to = "data.frame", def=function(from) {
    d <- .phylo4ToDataFrame(from, edgeOrder="pretty")
    d
})
