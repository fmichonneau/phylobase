## drop.tip4.R

## modified from drop.tip.R (2006-10-12)
##   Remove Tips in a Phylogenetic Tree
## Copyright 2003-2006 Emmanuel Paradis
## This file is part of the R-package `ape'.

## See the file ../COPYING for licensing issues.

setGeneric("prune",function(phy, ...) {
  standardGeneric("prune")
})

## setGeneric("drop.tip") ## if ape has already been loaded


DropTip <- function(phy,tip,...) {
  if (length(tip)==0) {
      phy
  } else if (is(phy,"phylo4d")) {
      ## use extract.tree instead of as() to avoid warning
      as(ape::drop.tip(as(extractTree(phy),"phylo"),tip,...),"phylo4d")
  } else as(ape::drop.tip(as(phy,"phylo"),tip,...),class(phy))
}


## return characters, sorted in NUMERIC order
.chnumsort <- function(x) {
  as.character(sort(as.numeric(x)))
}

setMethod("prune","phylo4",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...) {

    if (subtree) {
        #stop("subtree option is not currently supported for phylo4")
        # do this for now (at least to allow examples to pass check)
        warning("subtree option is not currently supported for phylo4")
        return(DropTip(phy,tip,trim.internal, subtree))
    }

    makeEdgeNames <- function(edge) {
        paste(edge[,1], edge[,2], sep="-")
    } 

    ## drop tips and obsolete internal nodes from edge matrix
    tip.drop <- getNode(phy, tip, missing="fail")
    tip.keep <- setdiff(seq_len(nTips(phy)), tip.drop)
    nodes <- seq_len(nrow(phy@edge))
    node.keep <- logical(nrow(phy@edge))
    node.keep[tip.keep] <- TRUE
    if (trim.internal) {
        if (phy@order == "postorder") {
            edge.post <- phy@edge
        } else {
            edge.post <- reorder(phy, "postorder")@edge
        }
        for (i in seq_along(edge.post[,2])) {
            if (node.keep[edge.post[i,2]]) {
                node.keep[edge.post[i,1]] <- TRUE
            }
        }
    } else {
        node.keep[nTips(phy) + seq_len(nNodes(phy))] <- TRUE
    }
    edge.new <- phy@edge[phy@edge[,2] %in% nodes[node.keep], ]
  
    ## remove singletons
    edge.length.new <- phy@edge.length
    edge.label.new <- phy@edge.label
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
    tip.label.new <- phy@tip.label[names(phy@tip.label) %in% edge.new]
    node.label.new <- phy@node.label[names(phy@node.label) %in% edge.new]

    ## subset and order edge.length and edge.label with respect to edge
    edge.names <- makeEdgeNames(edge.new)
    edge.length.new <- edge.length.new[edge.names]
    edge.label.new <- edge.label.new[edge.names]

    if (!trim.internal) {
        ## make sure now-terminal internal nodes are treated as tips
        tip.now <- setdiff(edge.new[,2], edge.new[,1])
        tip.add <- tip.now[tip.now>nTips(phy)]
        if (length(tip.add)>0) {
            ind <- match(tip.add, names(node.label.new))

            ## node renumbering workaround to satisfy plot method
            newid <- sapply(tip.add, function(x) descendants(phy, x)[1])
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
    ## NOTE: a faster but looser approach would be to replace the phy
    ## slots with their new values (including Nnode) and return phy
    phylo4(x=edge.new, edge.length = edge.length.new, tip.label =
        tip.label.new, node.label = node.label.new, edge.label =
        edge.label.new, annote=phy@annote)
})

## trace("prune", browser, signature = "phylo4d")
## untrace("prune", signature = "phylo4d")
setMethod("prune", "phylo4d", function(phy, tip, trim.internal=TRUE,
                                       subtree=FALSE, ...) {
    tree <- extractTree(phy)
    phytr <- prune(tree, tip, trim.internal, subtree)

    ## create temporary phylo4 object with unique labels
    tmpLbl <- .genlab("n", nTips(phy)+nNodes(phy))
    tmpPhy <- tree
    labels(tmpPhy, "all") <- tmpLbl
    tmpPhytr <- prune(tmpPhy, getNode(phy, tip), trim.internal, subtree)

    ## get node numbers to keep
    oldLbl <- labels(tmpPhy, "all")
    newLbl <- labels(tmpPhytr, "all")
    toKeep <- as.numeric(names(oldLbl[oldLbl %in% newLbl]))
    tipToKeep <- toKeep[toKeep %in% nodeId(phy, "tip")]
    nodToKeep <- toKeep[toKeep %in% nodeId(phy, "internal")]

    if(!all(dim(phy@tip.data) == 0)) {
        tipDt <- phy@tip.data[match(tipToKeep, rownames(phy@tip.data)) ,, drop=FALSE]
        tipDt <- tipDt[.chnumsort(rownames(tipDt)) ,, drop=FALSE]
        rownames(tipDt) <- 1:nTips(phytr)
    }
    else
        tipDt <- data.frame(NULL)

    if(!all(dim(phy@node.data) == 0)) {
        nodDt <- phy@node.data[match(nodToKeep, rownames(phy@node.data)) ,, drop=FALSE]
        nodDt <- nodDt[.chnumsort(rownames(nodDt)) ,, drop=FALSE]
        rownames(nodDt) <- 1:nNodes(phytr)
    }
    else
        nodDt <- data.frame(NULL)

    phytr <- phylo4d(phytr, tip.data=tipDt, node.data=nodDt, match.data=FALSE)

    phytr
})

setMethod("prune", "phylo",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   ...) {
            DropTip(phy, tip, trim.internal, subtree)
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

