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
             check1 <- check_data(object)
             check2 <- check_phylo4(object)          
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
setMethod("phylo4d", c("phylo4"), function(x, tip.data=NULL, node.data=NULL, all.data=NULL, ...){

    if(is.character(checkval <- check_phylo4(x))) stop(checkval)

    res <- new("phylo4d")
    res@edge <- x@edge
    res@edge.length <- x@edge.length
    res@Nnode <- x@Nnode
    res@tip.label <- x@tip.label
    res@node.label <- x@node.label
    res@edge.label <- x@edge.label
    res@root.edge <- x@root.edge

### handle a which argument
    which.dat <- match.arg(list(...)$"which", c("tip","node","all"))

    ## handle data
    if(all(is.null(c(tip.data, node.data, all.data)))) {
        stop("no data provided; please use phylo4 class")
    }

    ## convert vector to data.frames
    if(is.vector(tip.data)) tip.data <- as.data.frame(tip.data)
    if(is.vector(node.data)) node.data <- as.data.frame(node.data)
    if(is.vector(all.data)) all.data <- as.data.frame(all.data)

    if(!is.null(all.data)){
        if(!is.data.frame(all.data)) stop("all.data must be a data.frame")
        tip.data <- all.data[1:nTips(x) , , drop=FALSE]
        node.data <- all.data[-(1:nTips(x)) , , drop=FALSE]
    }

    ## now at least one data.frame is provided
    if(is.null(tip.data)) tip.data <- data.frame(NULL)
    if(is.null(node.data)) node.data <- data.frame(NULL)
    if(!is.data.frame(tip.data)) stop("tip.data must be a data.frame")
    if(!is.data.frame(node.data)) stop("node.data must be a data.frame")

    res@tip.data <- tip.data
    res@node.data <- node.data

    check_data(res, ...)
    res <- attach_data(res,...)
    return(res)  
})
