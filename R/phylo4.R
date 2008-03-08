require(methods)
require(ape)
         
## accessor functions for all internal bits
## HORRIBLE KLUGE
nTips <- function(x,...)  { }  ## mask ape::nTips
setGeneric("nTips", function(x,...) {
    standardGeneric("nTips")
})

## hack to ensure ape compatibility
setMethod("nTips","ANY", function(x) {
    if (class(x)=="phylo") {
        Ntip(x)
    } else stop(paste("no 'nTips' method available for",
                      deparse(substitute(x)),
                      "(class",class(x),")"))
})

setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})

setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})

setGeneric("edges", function(x,order,...) {
    standardGeneric("edges")
})

setGeneric("rootEdge", function(x,order,...) {
    standardGeneric("rootEdge")
})

setGeneric("isRooted", function(x) {
    standardGeneric("isRooted")
})

setGeneric("rootNode", function(x) {
    standardGeneric("rootNode")
})

setGeneric("rootNode<-", function(x,value) {
    standardGeneric("rootNode<-")
})

setGeneric("hasEdgeLength", function(x) {
    standardGeneric("hasEdgeLength")
})

setGeneric("edgeLength", function(x) {
    standardGeneric("edgeLength")
})

setGeneric("hasNodeLabels", function(x) {
    standardGeneric("hasNodeLabels")
})

setGeneric("hasEdgeLabels", function(x) {
    standardGeneric("hasEdgeLabels")
})

setGeneric("labels")

setGeneric("labels<-",
           function(object,...,value) {
               standardGeneric("labels<-")
           })

setGeneric("nodeLabels", function(x) {
    standardGeneric("nodeLabels")
})
setGeneric("nodeLabels<-",
           function(object,...,value) {
               standardGeneric("nodeLabels<-")
           })

setGeneric("edgeLabels", function(x) {
    standardGeneric("edgeLabels")
})

setGeneric("edgeLabels<-",
           function(object,...,value) {
               standardGeneric("edgeLabels<-")
           })

setGeneric("print")

setGeneric("tdata", function(x,...) {
    standardGeneric("tdata")
})

setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})

setGeneric("na.omit")

###################
## Function .genlab
###################
## recursive function to have labels of constant length
## base = a character string
## n = number of labels
.genlab <- function(base,n) {
    if (n<=0) return("")
    s <- seq(length.out=n)
    fw <- max(nchar(as.character(s)))
    numstr <- formatC(s,flag="0",width=fw)
    paste(base,numstr,sep="")
}

## convert from phylo to phylo4
## coerce phylo4d to phylo4 -- on purpose, so no warning
extract.tree <- function(from) {
    phylo4(edge=from@edge,
           edge.length=from@edge.length,
           Nnode=from@Nnode,
           tip.label=from@tip.label)
}

