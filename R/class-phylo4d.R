###################################
## phylo4d class
## extend: phylo with data
setClass("phylo4d",
         representation(data="data.frame",
                        metadata = "list"),

         prototype = list(
           data = data.frame(NULL),
           metadata = list()),

         validity = checkPhylo4,
         contains="phylo4")

######################
## phylo4d constructor
######################
## TEST ME
## '...' recognized args for data are tipdata and nodedata.
## other recognized options are those known by the phylo4 constructor
##

## generic
setGeneric("phylo4d", function(x, ...) { standardGeneric("phylo4d")} )

## Core part that takes care of the data
.phylo4Data <- function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                        merge.data=TRUE) {

    ## Check validity of phylo4 object
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)

    ## Create placeholder data frames for any null data arguments
    if (is.null(tip.data)) tip.data <- formatData(x, NULL, "tip")
    if (is.null(node.data)) node.data <- formatData(x, NULL, "internal")
    if (is.null(all.data)) all.data <- formatData(x, NULL, "all")

    # don't allow all.data columns of same name as tip.data or node.data
    colnamesTipOrNode <- union(names(tip.data), names(node.data))
    if (any(names(all.data) %in% colnamesTipOrNode)) {
        stop("all.data column names must be distinct from ",
             "tip.data and node.data column names")
    }

    ## combine common columns and move into all.data if merging,
    ## otherwise rename them
    colsToMerge <- intersect(names(tip.data), names(node.data))
    if (merge.data && length(colsToMerge)>0) {
        ##TODO could really just index rows directly on 1:nTip and
        ## (nTip+1):(nTip+nNode) in the next two statements for speed,
        ## but this is more robust to changes in node numbering rules
        tip.rows <- tip.data[match(nodeId(x, "tip"),
            row.names(tip.data)), colsToMerge, drop=FALSE]
        node.rows <- node.data[match(nodeId(x, "internal"),
            row.names(tip.data)), colsToMerge, drop=FALSE]
        merge.data <- rbind(tip.rows, node.rows)
        all.data <- data.frame(all.data, merge.data)
    } else {
        names(tip.data)[names(tip.data) %in% colsToMerge] <-
            paste(colsToMerge, "tip", sep=".")
        names(node.data)[names(node.data) %in% colsToMerge] <-
            paste(colsToMerge, "node", sep=".")
    }
    ## now separate tips-only and nodes-only data
    tip.only.data <- tip.data[setdiff(names(tip.data), names(node.data))]
    node.only.data <- node.data[setdiff(names(node.data), names(tip.data))]

    ## combine all data
    complete.data <- data.frame(all.data, tip.only.data, node.only.data)

    ## drop any rows that only contain NAs
    if (ncol(complete.data)==0) {
        return(data.frame())
    } else {
        empty.rows <- as.logical(rowSums(!is.na(complete.data)))
        return(complete.data[empty.rows, , drop=FALSE])
    }

}


## first arg is a phylo4
### phylo4d class rewrite
setMethod("phylo4d", "phylo4",
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                   merge.data=TRUE, metadata = list(), ...) {
    ## coerce tree to phylo4d
    res <- as(x, "phylo4d")

    ## apply formatData to ensure data have node number rownames and
    ## correct dimensions
    tip.data <- formatData(phy=x, dt=tip.data, type="tip", ...)
    node.data <- formatData(phy=x, dt=node.data, type="internal", ...)
    all.data <- formatData(phy=x, dt=all.data, type="all", ...)

    ## add any data
    res@data <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
        all.data=all.data, merge.data=merge.data)
    ## add any metadata
    res@metadata <- metadata
    return(res)
})


## first arg is a matrix of edges
setMethod("phylo4d", c("matrix"),
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                   merge.data=TRUE, metadata=list(), edge.length=NULL,
                   tip.label=NULL, node.label=NULL, edge.label=NULL,
                   order="unknown", annote=list(), ...) {
    tree <- phylo4(x, edge.length=edge.length, tip.label=tip.label,
        node.label=node.label, edge.label=edge.label, order=order,
        annote=annote)
    res <- phylo4d(tree, tip.data, node.data, all.data,
        merge.data=merge.data, metadata=metadata, ...)
    return(res)
})

## first arg is a phylo
setMethod("phylo4d", "phylo",
          function(x, tip.data=NULL,
                   node.data=NULL, all.data=NULL,
                   check.node.labels=c("keep", "drop", "asdata"),
                   annote=list(), metadata=list(), ...) {

    check.node.labels <- match.arg(check.node.labels)

    if (check.node.labels == "asdata") {
        # FIXME? use.node.names=TRUE won't work with this option b/c
        # node labels are dropped; assumes node.data (if any), phylo
        # node.label, and phylo4 internal nodes are in the same order?

        nlab.data <- x$node.label
        x$node.label <- NULL
        nlab.data[!nzchar(nlab.data)] <- NA

        nlab.data <- data.frame(labelValues=as.numeric(nlab.data))

        tree <- phylo4(x, check.node.labels="drop", annote=annote)
        res <- phylo4d(tree, tip.data=tip.data, node.data=node.data,
                       all.data=all.data, metadata=metadata, ...)
        res <- addData(res, node.data=nlab.data, pos="before", match.data=FALSE)
    }
    else {
        tree <- phylo4(x, check.node.labels=check.node.labels, annote=annote)
        res <- phylo4d(tree, tip.data=tip.data, node.data=node.data,
                       all.data=all.data, metadata=metadata, ...)
    }

    return(res)
})

### first arg is a phylo4d
setMethod("phylo4d", c("phylo4d"), function(x, ...) {
          stop("Your object is already a phylo4d object. If you want to modify the data attached to it look at the help for tdata()<-")
      })

