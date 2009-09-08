## Same order as in methods-phylo4.R

## nTips
setGeneric("nTips", function(x,...) {
    standardGeneric("nTips")
})

## nNodes
setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})

## nodeType
setGeneric("nodeType", function(phy) {
    standardGeneric("nodeType")
})

## nodeId
setGeneric("nodeId", function(x, type=c("all", "tip", "internal",
    "root")) {
    standardGeneric("nodeId")
})

## nEdges
setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})

## edges
setGeneric("edges", function(x, order, ...) {
    standardGeneric("edges")
})

## edgeOrder
setGeneric("edgeOrder", function(x, ...) {
  standardGeneric("edgeOrder")
})

## edgeId
setGeneric("edgeId", function(x, type=c("all", "tip", "internal",
    "root")) {
    standardGeneric("edgeId")
})

## hasEdgeLength
setGeneric("hasEdgeLength", function(x) {
    standardGeneric("hasEdgeLength")
})

## edgeLength
setGeneric("edgeLength", function(x, ...) {
    standardGeneric("edgeLength")
})

## edgeLength<-
setGeneric("edgeLength<-", function(x, ..., value) {
    standardGeneric("edgeLength<-")
})

## sumEdgeLength
setGeneric("sumEdgeLength", function(phy, node) {
    standardGeneric("sumEdgeLength")
})

## isRooted
setGeneric("isRooted", function(x) {
    standardGeneric("isRooted")
})

## rootNode
setGeneric("rootNode", function(x) {
    standardGeneric("rootNode")
})

## rootNode<-
setGeneric("rootNode<-", function(x, value) {
    standardGeneric("rootNode<-")
})

## labels
setGeneric("labels")

## labels<-
setGeneric("labels<-",
           function(object, type, use.names, ..., value) {
               standardGeneric("labels<-")
           })

## hasNodeLabels
setGeneric("hasNodeLabels", function(x) {
    standardGeneric("hasNodeLabels")
})

## nodeLabels
setGeneric("nodeLabels", function(object) {
    standardGeneric("nodeLabels")
})

## nodeLabels<-
setGeneric("nodeLabels<-",
           function(object, ..., value) {
               standardGeneric("nodeLabels<-")
           })

## tipLabels
setGeneric("tipLabels", function(object) {
    standardGeneric("tipLabels")
})

## tipLabels<-
setGeneric("tipLabels<-",
   function(object, ..., value) {
       standardGeneric("tipLabels<-")
   })

## hasEdgeLabels
setGeneric("hasEdgeLabels", function(x) {
    standardGeneric("hasEdgeLabels")
})

## edgeLabels
setGeneric("edgeLabels", function(x) {
    standardGeneric("edgeLabels")
})

## edgeLabels<-
setGeneric("edgeLabels<-",
           function(object, ..., value) {
               standardGeneric("edgeLabels<-")
           })

## print
setGeneric("print")

## head
setGeneric("head")

## tail
setGeneric("tail")

### ----------- phylo4d methods -----------

## tdata
setGeneric("tdata", function(x, ...) {
    standardGeneric("tdata")
})

## tdata<-
setGeneric("tdata<-", function(x, ..., value) {
    standardGeneric("tdata<-")
})

## addData
setGeneric("addData", function(x, ...) {
    standardGeneric("addData")
})

## hasTipData
setGeneric("hasTipData", function(x) {
    standardGeneric("hasTipData")
})

## hasNodeData
setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})

##
##setGeneric("na.omit")

setGeneric("reorder")

##setGeneric("rootEdge", function(x, order, ...) {
##    standardGeneric("rootEdge")
##})

###################
## Function .genlab
###################
## (formerly) recursive function to have labels of constant length
## base = a character string
## n = number of labels
.genlab <- function(base, n) {
    if(n <= 0) return("")
    s <- seq(length.out=n)
    fw <- max(nchar(as.character(s)))
    numstr <- formatC(s, flag="0", width=fw)
    paste(base, numstr, sep="")
}


