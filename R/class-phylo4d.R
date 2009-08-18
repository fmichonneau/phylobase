###################################
## phylo4d class
## extend: phylo with data
setClass("phylo4d",
         representation(tip.data="data.frame",
                        node.data="data.frame"),
         ##                        edgedata="data.frame"),
         prototype = list( tip.data = data.frame(NULL),
           node.data = data.frame(NULL) ),
         ##all.data = data.frame(NULL) ),
         validity = function(object) {
             ## FIXME: finish this by intercepting FALSE, char string, etc.
             check1 <- checkData(object)
             check2 <- checkPhylo4(object)
         },
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
                        match.data=TRUE, merge.data=TRUE, ...) {

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

    ## Replacing node labels by node numbers and formatting the data to make sure
    ## they have the correct dimensions
    if(!is.null(all.data) && all(dim(all.data) > 0))
        all.data <- formatData(x, all.data, which="all",
                               match.data=match.data, ...)

    if(!is.null(tip.data) && all(dim(tip.data) > 0))
        tip.data <- formatData(x, tip.data, which="tip",
                               match.data=match.data, ...)

    if(!is.null(node.data) && all(dim(node.data) > 0))
        node.data <- formatData(x, node.data, which="internal",
                                match.data=match.data, ...)

    ## Merging dataset
    if(!is.null(all.data)) {
        tmpData <- all.data
        if(!is.null(tip.data)) {
            emptyNodeData <- array(, dim = c(nNodes(x), ncol(tip.data)),
                                   dimnames = list(nodeId(x, "internal"),
                                   colnames(tip.data)))
            tmpTipData <- rbind(tip.data, emptyNodeData)

            tmpTipData <- tmpTipData[match(rownames(all.data),
                                           rownames(tmpTipData)) ,,
                                     drop = FALSE]
            tmpData <- cbind(all.data, tmpTipData)
        }
        if(!is.null(node.data)) {
            emptyTipData <- array(, dim = c(nTips(x), ncol(node.data)),
                                  dimnames = list(nodeId(x, "tip"),
                                  colnames(node.data)))
            tmpNodeData <- rbind(emptyTipData, node.data)
            tmpNodeData <- tmpNodeData[match(rownames(all.data),
                                             rownames(tmpNodeData)) ,,
                                       drop = FALSE]
            tmpData <- cbind(tmpData, tmpNodeData)
        }

        if(match.data) {
            tip.data <- tmpData[rownames(tmpData) %in% nodeId(x, "tip") ,,
                                    drop = FALSE]
            node.data <- tmpData[rownames(tmpData) %in% nodeId(x, "internal") ,,
                                     drop = FALSE]
        }
        else {
            tip.data <- tmpData[1:nTips(x) ,, drop=FALSE]
            node.data <- tmpData[-(1:nTips(x)) ,, drop=FALSE]
        }

    }

    else {
        if(!is.null(tip.data) && !is.null(node.data)) {
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
            if(is.null(tip.data)) tip.data <- data.frame(NULL)
            if(is.null(node.data)) node.data <- data.frame(NULL)

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
                   match.data=TRUE, merge.data=TRUE, ...) {

    ## Creating new phylo4d object
    res <- new("phylo4d")
    res@edge <- x@edge
    res@edge.length <- x@edge.length
    res@Nnode <- x@Nnode
    res@tip.label <- x@tip.label
    res@node.label <- x@node.label
    res@edge.label <- x@edge.label

    ## taking care of the data
    tmpData <- .phylo4Data(x, tip.data, node.data, all.data, match.data,
                           merge.data, ...)

    res@tip.data <- tmpData$tip.data
    res@node.data <- tmpData$node.data

    return(res)
})


## first arg is a matrix of edges
setMethod("phylo4d", c("matrix"),
          function(x, tip.data=NULL, node.data=NULL, all.data=NULL, ...) {
    tree <- phylo4(x, ...)
    res <- phylo4d(tree, tip.data, node.data, all.data, ...)
    return(res)
})

## first arg is a phylo
setMethod("phylo4d", c("phylo"),
          function(x, tip.data=NULL,
                   node.data=NULL, all.data=NULL,
                   check.node.labels=c("keep", "drop", "asdata"), ...) {

    check.node.labels <- match.arg(check.node.labels)

    if (check.node.labels == "asdata") {
        # FIXME? use.node.names=TRUE won't work with this option b/c
        # node labels are dropped; assumes node.data (if any), phylo
        # node.label, and phylo4 internal nodes are in the same order?

        nlab.data <- x$node.label
        x$node.label <- NULL
        nlab.data[!nzchar(nlab.data)] <- NA

        nlab.data <- data.frame(labelValues=as.numeric(nlab.data))

        tree <- phylo4(x, check.node.labels="drop")
        res <- phylo4d(tree, tip.data, node.data, all.data, ...)
        res <- addData(res, node.data=nlab.data, pos="before", match.data=FALSE)
    }
    else {
        tree <- phylo4(x, check.node.labels=check.node.labels)
        res <- phylo4d(tree, tip.data, node.data, all.data, ...)
    }

    return(res)
})

### first arg is a phylo4d
setMethod("phylo4d", c("phylo4d"), function(x, ...) {
          stop("Your object is already a phylo4d object. If you want to modify the data attached to it look at the help for tdata()<-")
      })

