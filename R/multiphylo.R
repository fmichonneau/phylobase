require(methods)
require(ape)

setClass("multiPhylo4",
         representation(phylolist="list",
                        tree.names="character"),
         prototype = list(phylolist=list(),
           tree.names=character(0)))

setClass("multiPhylo4d",
         representation(tip.data="data.frame"),
         contains="multiPhylo4")


## function to bind trees together into a multi-tree object
tbind <- function(...,check_data=FALSE) {
    L <- as.list(...)
    treeclasses <- c("multiPhylo4d","multiPhylo4","phylo4","phylo4d")
    tdataclasses <- c("multiPhylo4d","phylo4d")
    classes <- sapply(L,class)
    if (!all(classes %in% treeclasses)) stop("all elements must be trees or multitrees")
    if (!all(classes %in% tdataclasses)) {
        if (any(classes %in% tdataclasses)) warning("not all elements contain data: data discarded")
        ## decompose multi-trees into lists
        ## bind list into multi-tree
    } else {
        ## check: all data identical?
        ## decompose multi-trees into lists
    }
}
