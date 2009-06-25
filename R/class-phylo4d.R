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

## first arg is a phylo4
setMethod("phylo4d", "phylo4",
   function(x, tip.data = NULL, node.data = NULL, all.data = NULL,
            merge.tip.node = TRUE, ...) {

       classData <- function(someData) {
           if(!is.null(someData)) {
               if(is.vector(someData)) someData <- as.data.frame(someData)
               if(!is.data.frame(someData)) {
                   nmSomedata <- deparseSubstitute(someData)
                   return(paste(nmSomeData, "must be a vector or a data frame"))
               }
               return(TRUE)
           }
           else return(TRUE)
       }

       if(is.character(checkval <- checkPhylo4(x))) stop(checkval)

       if(is.character(checkClass <- classData(all.data))) stop(checkClass)
       if(is.character(checkClass <- classData(tip.data))) stop(checkClass)
       if(is.character(checkClass <- classData(node.data))) stop(checkClass)

       res <- new("phylo4d")
       res@edge <- x@edge
       res@edge.length <- x@edge.length
       res@Nnode <- x@Nnode
       res@tip.label <- x@tip.label
       res@node.label <- x@node.label
       res@edge.label <- x@edge.label

       if(!is.null(all.data)) {
           tmpData <- all.data
           if(!is.null(tip.data)) {
               emptyNodeData <- array(, dim = c(nNodes(x), ncol(tip.data)),
                                      dimnames = list(nodeLabels(x), colnames(tip.data)))
               tmpTipData <- rbind(tip.data, emptyNodeData)
               ## TODO? - have a test on names between
               tmpTipData <- tmpTipData[match(rownames(all.data), rownames(tmpTipData)) ,, drop = FALSE]
               tmpData <- cbind(all.data, tmpTipData)
           }
           if(!is.null(node.data)) {
               emptyTipData <- array(, dim = c(nTips(x), ncol(node.data)),
                                     dimnames = list(tipLabels(x), colnames(node.data)))
               tmpNodeData <- rbind(emptyTipData, node.data)
               ## TODO? - add test
               tmpNodeData <- tmpNodeData[match(rownames(all.data), rownames(tmpNodeData)) ,, drop = FALSE]
               tmpData <- cbind(tmpData, tmpNodeData)

           }
           if (!hasNodeLabels(x)) stop("can't match node data to labels without node labels")
           res@tip.data <- tmpData[rownames(tmpData) %in% tipLabels(x) ,, drop = FALSE]
           res@node.data <- tmpData[rownames(tmpData) %in% nodeLabels(x) ,, drop = FALSE]
       }

       else {
           if((!is.null(tip.data) && (!is.null(node.data)))) {
               if(identical(colnames(tip.data), colnames(node.data)) &&  merge.tip.node) {
                   tmpAllData <- rbind(tip.data, node.data)
                   res@tip.data <- tmpAllData[1:nTips(x) ,, drop = FALSE]
                   res@node.data <- tmpAllData[-(1:nTips(x)) ,, drop = FALSE]
               }
               else {
                   emptyTipData <- array(, dim = c(nTips(x), ncol(node.data)),
                                           dimnames = list(tipLabels(x), colnames(node.data)))
                   emptyNodeData <- array(, dim = c(nNodes(x), ncol(tip.data)),
                                            dimnames = list(nodeLabels(x), colnames(tip.data)))
                   tmpTipData <- rbind(tip.data, emptyNodeData)
                   tmpNodeData <- rbind(emptyTipData, node.data)
                   tmpData <- cbind(tmpTipData, tmpNodeData)
                   res@tip.data <- tmpData[1:nTips(x) ,, drop = FALSE]
                   res@node.data <- tmpData[-(1:nTips(x)) ,, drop = FALSE]
               }
           }
           else {
               ## at this point provide NULL data frame for empty arguments
               if(is.null(tip.data)) tip.data <- data.frame(NULL)
               if(is.null(node.data)) node.data <- data.frame(NULL)

               res@tip.data <- tip.data
               res@node.data <- node.data
           }
       }

       checkData(res, ...)
       res <- attachData(res,...)
       return(res)

})

## first arg is a matrix of edges
setMethod("phylo4d", c("matrix"), function(x, tip.data=NULL, node.data=NULL, all.data=NULL, ...){
    tree <- phylo4(x, ...)
    res <- phylo4d(tree, tip.data, node.data, all.data, ...)
    return(res)
})

## first arg is a phylo
setMethod("phylo4d", c("phylo"), function(x, tip.data=NULL,
    node.data=NULL, all.data=NULL, check.node.labels=c("keep", "drop",
    "asdata"), ...) {

    check.node.labels <- match.arg(check.node.labels)
    if (check.node.labels == "asdata") {
        # FIXME? use.node.names=TRUE won't work with this option b/c
        # node labels are dropped; assumes node.data (if any), phylo
        # node.label, and phylo4 internal nodes are in the same order?
        nlab.data <- x$node.label
        x$node.label <- NULL
        nlab.data[!nzchar(nlab.data)] <- NA
        # TODO only convert to numeric if values are number-like?
        nlab.data <- data.frame(labelValues=as.numeric(nlab.data))
        if (is.null(node.data)) {
            node.data <- nlab.data
        } else {
            node.data <- cbind(nlab.data, node.data)
        }
        tree <- phylo4(x, check.node.labels="drop")
    } else {
        tree <- phylo4(x, check.node.labels=check.node.labels)
    }
    res <- phylo4d(tree, tip.data, node.data, all.data, ...)
    return(res)
})

### first arg is a phylo4d
setMethod("phylo4d", c("phylo4d"), function(x, ...) {
          stop("Your object is already a phylo4d object. If you want to modify the data attached to it look at the help for tdata()<-")
      })

