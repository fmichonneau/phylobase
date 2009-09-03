###################################
## phylo4d class
## extend: phylo with data
setClass("phylo4d",
         representation(tip.data="data.frame",
                        node.data="data.frame",
                        metadata = "list"),

         prototype = list( tip.data = data.frame(NULL),
           node.data = data.frame(NULL),
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
                        match.data=TRUE, merge.data=TRUE,
                        rownamesAsLabels=FALSE,
                        ...) {

    ## Make sure that data provided are a data frame
    classData <- function(someData) {
        if(!is.null(someData)) {
            if(is.vector(someData))
                someData <- as.data.frame(someData)
            if(!is.data.frame(someData)) {
                nmSomeData <- substitute(someData)
                stop(paste(nmSomeData, "must be a vector or a data frame"))
            }
            someData
        }
    }

    ## Check validity of phylo4 object
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)

    ## Check/Transform provided data to data.frame
    all.data <- classData(all.data)
    tip.data <- classData(tip.data)
    node.data <- classData(node.data)

    is.empty <- function(x) { is.null(x) || all(dim(x)==0) }

    ## Replacing node labels by node numbers and formatting the data to make sure
    ## they have the correct dimensions
    if(!is.empty(all.data))
        all.data <- formatData(phy=x, dt=all.data, type="all",
                               match.data=match.data,
                               rownamesAsLabels=rownamesAsLabels, ...)

    if(!is.empty(tip.data))
        tip.data <- formatData(phy=x, dt=tip.data, type="tip",
                               match.data=match.data,
                               rownamesAsLabels=rownamesAsLabels, ...)

    if(!is.empty(node.data))
        node.data <- formatData(phy=x, dt=node.data, type="internal",
                                match.data=match.data,
                                rownamesAsLabels=rownamesAsLabels, ...)

    ## Merging dataset
    if(!is.empty(all.data)) {
        tmpData <- all.data
        if(!is.empty(tip.data)) {
            emptyNodeData <- array(, dim = c(nNodes(x), ncol(tip.data)),
                                   dimnames = list(nodeId(x, "internal"),
                                   colnames(tip.data)))
            tmpTipData <- rbind(tip.data, emptyNodeData)

            tmpTipData <- tmpTipData[match(rownames(all.data),
                                           rownames(tmpTipData)) ,,
                                     drop = FALSE]
            tmpData <- cbind(all.data, tmpTipData)
        }
        if(!is.empty(node.data)) {
            emptyTipData <- array(, dim = c(nTips(x), ncol(node.data)),
                                  dimnames = list(nodeId(x, "tip"),
                                  colnames(node.data)))
            tmpNodeData <- rbind(emptyTipData, node.data)
            tmpNodeData <- tmpNodeData[match(rownames(all.data),
                                             rownames(tmpNodeData)) ,,
                                       drop = FALSE]
            tmpData <- cbind(tmpData, tmpNodeData)
        }

        tip.data <- tmpData[rownames(tmpData) %in% nodeId(x, "tip") ,,
                            drop = FALSE]
        node.data <- tmpData[rownames(tmpData) %in% nodeId(x, "internal") ,,
                             drop = FALSE]
    }

    else {
        if(!is.empty(tip.data) && !is.empty(node.data)) {
            if(identical(colnames(tip.data), colnames(node.data)) && merge.data) {
                tmpAllData <- rbind(tip.data, node.data)
                tip.data <- tmpAllData[rownames(tmpAllData) %in%
                                           nodeId(x, "tip") ,, drop=FALSE]
                node.data <- tmpAllData[rownames(tmpAllData) %in%
                                            nodeId(x, "internal") ,, drop=FALSE]
            }
            else {
                emptyTipData <- array(, dim = c(nTips(x), ncol(node.data)),
                                      dimnames = list(nodeId(x, "tip"),
                                      colnames(node.data)))
                emptyNodeData <- array(, dim = c(nNodes(x), ncol(tip.data)),
                                       dimnames = list(nodeId(x, "internal"),
                                       colnames(tip.data)))
                tmpTipData <- rbind(tip.data, emptyNodeData)
                tmpNodeData <- rbind(emptyTipData, node.data)
                tmpNodeData <- tmpNodeData[rownames(tmpTipData) ,, drop=FALSE]

                tmpData <- cbind(tmpTipData, tmpNodeData)

                if(match.data) {
                    tip.data <- tmpData[rownames(tmpData) %in%
                                            nodeId(x, "tip") ,, drop=FALSE]
                    node.data <- tmpData[rownames(tmpData) %in%
                                             nodeId(x, "internal") ,, drop=FALSE]
                }
                else {
                    tip.data <- tmpData[1:nTips(x) ,, drop=FALSE]
                    node.data <- tmpData[-(1:nTips(x)) ,, drop=FALSE]
                }
            }
        }
        else {
            ## at this point provide NULL data frame for empty arguments
            if(is.empty(tip.data)) tip.data <- data.frame(NULL)
            if(is.empty(node.data)) node.data <- data.frame(NULL)

            tip.data <- tip.data
            node.data <- node.data
        }
    }

    return(list(tip.data=tip.data, node.data=node.data))
}


## first arg is a phylo4
### phylo4d class rewrite
setMethod("phylo4d", "phylo4",
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                   match.data=TRUE, merge.data=TRUE, rownamesAsLabels=FALSE,
                   metadata = list(),
                   ...) {

    ## prepare the data
    tmpData <- .phylo4Data(x=x, tip.data=tip.data, node.data=node.data,
                           all.data=all.data, match.data=match.data,
                           merge.data=merge.data,
                           rownamesAsLabels=rownamesAsLabels, ...)

    ## coerce to phylo4d and add data/metadata
    res <- as(x, "phylo4d")
    res@tip.data <- tmpData$tip.data
    res@node.data <- tmpData$node.data
    res@metadata <- metadata
    return(res)
})


## first arg is a matrix of edges
setMethod("phylo4d", c("matrix"),
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
                   metadata = list(), ...) {
    tree <- phylo4(x, ...)
    res <- phylo4d(tree, tip.data, node.data, all.data, metadata, ...)
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

