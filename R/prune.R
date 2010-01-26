
setGeneric("prune", function(x, ...) {
    standardGeneric("prune")
})


## return characters, sorted in NUMERIC order
.chnumsort <- function(x) {
  as.character(sort(as.numeric(x)))
}

setMethod("prune","phylo4", function(x, tips.exclude,
    trim.internal=TRUE) {

    makeEdgeNames <- function(edge) {
        paste(edge[,1], edge[,2], sep="-")
    } 

    ## drop tips and obsolete internal nodes from edge matrix
    tip.drop <- getNode(x, tips.exclude, missing="fail")
    tip.keep <- setdiff(nodeId(x, "tip"), tip.drop)
    nodes <- nodeId(x, "all")
    node.keep <- rep(FALSE, length(nodes))
    node.keep[tip.keep] <- TRUE
    if (trim.internal) {
        if (edgeOrder(x) == "postorder") {
            edge.post <- edges(x)
        } else {
            browser()
            edge.post <- edges(reorder(x, "postorder"))
        }
        for (i in seq_along(edge.post[,2])) {
            if (node.keep[edge.post[i,2]]) {
                node.keep[edge.post[i,1]] <- TRUE
            }
        }
    } else {
        node.keep[nodeId(x, "internal")] <- TRUE
    }
    edge.new <- edges(x)[edges(x)[,2] %in% nodes[node.keep], ]
  
    ## remove singletons
    edge.length.new <- edgeLength(x)
    edge.label.new <- edgeLabels(x)
    singletons <- which(tabulate(edge.new[edge.new[, 1] != 0, 1])==1)
    while (length(singletons)>0) {
        sing.node <- singletons[1]

        ## update edge matrix
        edges.drop <- which(edge.new==sing.node, arr.ind=TRUE)[,"row"]
        sing.edges <- edge.new[edges.drop,]
        edge.new[edges.drop[2], ] <- c(sing.edges[2,1], sing.edges[1,2])
        edge.new <- edge.new[-edges.drop[1], ]

        ## update edge lengths and edge labels
        edge.names.drop <- makeEdgeNames(sing.edges)
        edge.name.new <- paste(sing.edges[2,1], sing.edges[1,2], sep="-")
        edge.length.new[edge.name.new] <-
            sum(edge.length.new[edge.names.drop])
        edge.length.new <- edge.length.new[-match(edge.names.drop,
            names(edge.length.new))]
        edge.label.new[edge.name.new] <- NA
        browser()
        edge.label.new <- edge.label.new[-match(edge.names.drop,
            names(edge.label.new))]

        singletons <- which(tabulate(edge.new[edge.new[, 1] != 0, 1])==1)
    }

    ## remove dropped elements from labels
    label.new <- labels(x)[names(labels(x)) %in% edge.new]

    ## subset and order edge.length and edge.label with respect to edge
    edge.names <- makeEdgeNames(edge.new)
    edge.length.new <- edge.length.new[edge.names]
    edge.label.new <- edge.label.new[edge.names]

    if (!trim.internal) {
        ## make sure now-terminal internal nodes are treated as tips
        tip.now <- setdiff(edge.new[,2], edge.new[,1])
        tip.add <- tip.now[tip.now>nTips(x)]
        if (length(tip.add)>0) {
            ind <- match(tip.add, names(label.new))

            ## node renumbering workaround to satisfy plot method
            newid <- sapply(tip.add, function(tip) descendants(x, tip)[1])
            names(label.new)[ind] <- newid
            edge.new[match(tip.add, edge.new)] <- newid
            tip.now[match(tip.add, tip.now)] <- newid

            isTip <- edge.new %in% tip.now
            edge.new[isTip] <- match(edge.new[isTip],
            sort(unique.default(edge.new[isTip])))
        }
    }

    ## renumber nodes in the edge matrix
    edge.new[] <- match(edge.new, sort(unique.default(edge.new))) - 1L

    ## update corresponding element names in the other slots
    edge.names <- makeEdgeNames(edge.new)
    names(edge.length.new) <- edge.names
    names(edge.label.new) <- edge.names
    label.new <- label.new[order(as.numeric(names(label.new)))]
    names(label.new) <- seq_along(label.new)

    ## update, check, then return the pruned phylo4 object
    x@edge <- edge.new
    ##TODO would prefer to leave out NA from edge.length slot, but can't 
    x@edge.length <- edge.length.new
    x@edge.label <- edge.label.new[!is.na(edge.label.new)]
    x@label <- label.new[!is.na(label.new)]
    if(is.character(checkval <- checkPhylo4(x))) {
        stop(checkval)
    } else {
        return(x)
    }

})

## trace("prune", browser, signature = "phylo4d")
## untrace("prune", signature = "phylo4d")
setMethod("prune", "phylo4d", function(x, tips.exclude,
    trim.internal=TRUE) {

    tree <- extractTree(x)
    phytr <- prune(tree, tips.exclude, trim.internal)

    ## create temporary phylo4 object with complete and unique labels
    tmpLbl <- .genlab("n", nTips(x)+nNodes(x))
    tmpPhy <- tree
    labels(tmpPhy, "all") <- tmpLbl
    tmpPhytr <- prune(tmpPhy, getNode(x, tips.exclude), trim.internal)

    ## get node numbers to keep
    oldLbl <- labels(tmpPhy, "all")
    newLbl <- labels(tmpPhytr, "all")
    wasKept <- oldLbl %in% newLbl
    nodesToKeep <- as.numeric(names(oldLbl[wasKept]))

    ## subset original data, and update names
    allDt <- x@data[match(nodesToKeep, rownames(x@data)), , drop=FALSE]
    rownames(allDt) <- match(newLbl, oldLbl[wasKept])

    phytr <- phylo4d(phytr, all.data=allDt, match.data=TRUE)

    phytr
})

## setMethod("prune","ANY",
##           function(phy, tip, trim.internal = TRUE, subtree = FALSE,
##                    ,...) {
##             if (class(phy)=="phylo") {
##               ape::prune(phy, tip, trim.internal, subtree)
##               } else stop("no prune method available for",
##                     deparse(substitute(phy)),
##                     "(class",class(phy),")")
##           })

