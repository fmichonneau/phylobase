
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
    singletons <- which(tabulate(na.omit(edge.new[,1]))==1)
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
        edge.label.new <- edge.label.new[-match(edge.names.drop,
            names(edge.label.new))]

        singletons <- which(tabulate(na.omit(edge.new[,1]))==1)
    }

    ## remove dropped elements from tip.label and node.label
    tip.label.new <- tipLabels(x)[names(tipLabels(x)) %in% edge.new]
    node.label.new <- nodeLabels(x)[names(nodeLabels(x)) %in% edge.new]

    ## subset and order edge.length and edge.label with respect to edge
    edge.names <- makeEdgeNames(edge.new)
    edge.length.new <- edge.length.new[edge.names]
    edge.label.new <- edge.label.new[edge.names]

    if (!trim.internal) {
        ## make sure now-terminal internal nodes are treated as tips
        tip.now <- setdiff(edge.new[,2], edge.new[,1])
        tip.add <- tip.now[tip.now>nTips(x)]
        if (length(tip.add)>0) {
            ind <- match(tip.add, names(node.label.new))

            ## node renumbering workaround to satisfy plot method
            newid <- sapply(tip.add, function(tip) descendants(x, tip)[1])
            names(node.label.new)[ind] <- newid
            edge.new[match(tip.add, edge.new)] <- newid
            tip.now[match(tip.add, tip.now)] <- newid

            tip.label.new <- c(tip.label.new, node.label.new[ind])
            node.label.new <- node.label.new[-ind]
            isTip <- edge.new %in% tip.now
            edge.new[isTip] <- match(edge.new[isTip],
            sort(unique.default(edge.new[isTip])))
        }
    }

    ## renumber nodes in the edge matrix
    edge.new[] <- match(edge.new, sort(unique.default(edge.new)))

    ## update corresponding element names in the other slots
    edge.names <- makeEdgeNames(edge.new)
    names(edge.length.new) <- edge.names
    names(edge.label.new) <- edge.names
    tip.label.new <- tip.label.new[order(as.numeric(names(tip.label.new)))]
    names(tip.label.new) <- seq_along(tip.label.new)
    names(node.label.new) <- seq_along(node.label.new) + length(tip.label.new)

    ## create and return new phylo4 object
    ## NOTE: a faster but looser approach would be to replace the slots
    ## of x with their new values (including Nnode) and return x
    phylo4(x=edge.new, edge.length = edge.length.new, tip.label =
        tip.label.new, node.label = node.label.new, edge.label =
        edge.label.new, annote=x@annote)
})

## trace("prune", browser, signature = "phylo4d")
## untrace("prune", signature = "phylo4d")
setMethod("prune", "phylo4d", function(x, tips.exclude,
    trim.internal=TRUE) {

    tree <- extractTree(x)
    phytr <- prune(tree, tips.exclude, trim.internal)

    ## create temporary phylo4 object with unique labels
    tmpLbl <- .genlab("n", nTips(x)+nNodes(x))
    tmpPhy <- tree
    labels(tmpPhy, "all") <- tmpLbl
    tmpPhytr <- prune(tmpPhy, getNode(x, tips.exclude), trim.internal)

    ## get node numbers to keep
    oldLbl <- labels(tmpPhy, "all")
    newLbl <- labels(tmpPhytr, "all")
    toKeep <- as.numeric(names(oldLbl[oldLbl %in% newLbl]))
    tipToKeep <- toKeep[toKeep %in% nodeId(x, "tip")]
    nodToKeep <- toKeep[toKeep %in% nodeId(x, "internal")]

    if(!all(dim(x@tip.data) == 0)) {
        tipDt <- x@tip.data[match(tipToKeep, rownames(x@tip.data)) ,, drop=FALSE]
        tipDt <- tipDt[.chnumsort(rownames(tipDt)) ,, drop=FALSE]
        rownames(tipDt) <- 1:nTips(phytr)
    }
    else
        tipDt <- data.frame(NULL)

    if(!all(dim(x@node.data) == 0)) {
        nodDt <- x@node.data[match(nodToKeep, rownames(x@node.data)) ,, drop=FALSE]
        nodDt <- nodDt[.chnumsort(rownames(nodDt)) ,, drop=FALSE]
        rownames(nodDt) <- 1:nNodes(phytr)
    }
    else
        nodDt <- data.frame(NULL)

    phytr <- phylo4d(phytr, tip.data=tipDt, node.data=nodDt, match.data=FALSE)

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

