## Same order as in methods-phylo4.R

## nTips
setGeneric("nTips", function(x) {
    standardGeneric("nTips")
})

## depthTips
setGeneric("depthTips", function(x) {
  standardGeneric("depthTips")
})

## nNodes
setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})

## nodeType
setGeneric("nodeType", function(x) {
    standardGeneric("nodeType")
})

## nodeId
setGeneric("nodeId", function(x, type=c("all", "tip", "internal",
    "root")) {
    standardGeneric("nodeId")
})

## nodeDepth
setGeneric("nodeDepth", function(x, node) {
  standardGeneric("nodeDepth")
})

## nEdges
setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})

## edges
setGeneric("edges", function(x, ...) {
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
setGeneric("sumEdgeLength", function(x, node) {
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
           function(x, type, use.names, ..., value) {
               standardGeneric("labels<-")
           })

## hasDuplicatedLabels
setGeneric("hasDuplicatedLabels",
           function(x, type) {
               standardGeneric("hasDuplicatedLabels")
           })

## hasNodeLabels
setGeneric("hasNodeLabels", function(x) {
    standardGeneric("hasNodeLabels")
})

## nodeLabels
setGeneric("nodeLabels", function(x) {
    standardGeneric("nodeLabels")
})

## nodeLabels<-
setGeneric("nodeLabels<-",
           function(x, ..., value) {
               standardGeneric("nodeLabels<-")
           })

## tipLabels
setGeneric("tipLabels", function(x) {
    standardGeneric("tipLabels")
})

## tipLabels<-
setGeneric("tipLabels<-",
   function(x, ..., value) {
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
           function(x, ..., value) {
               standardGeneric("edgeLabels<-")
           })

## print
setGeneric("print")


## head
setGeneric("head")

## tail
setGeneric("tail")

## summary
setGeneric("summary")

## isUltrametric
setGeneric("isUltrametric", function(x, tol=.Machine$double.eps^.5) {
  standardGeneric("isUltrametric")
})


### ----------- phylo4d methods -----------

## tdata
setGeneric("tdata", function(x, ...) {
    standardGeneric("tdata")
})

## tdata<-
setGeneric("tdata<-", function(x, ..., value) {
    standardGeneric("tdata<-")
})

## tipData
setGeneric("tipData", function(x, ...) {
    standardGeneric("tipData")
})

## tipData<-
setGeneric("tipData<-", function(x, ..., value) {
    standardGeneric("tipData<-")
})

## nodeData
setGeneric("nodeData", function(x, ...) {
    standardGeneric("nodeData")
})

## nodeData<-
setGeneric("nodeData<-", function(x, ..., value) {
    standardGeneric("nodeData<-")
})

## addData
setGeneric("addData", function(x, ...) {
    standardGeneric("addData")
})

## nData
setGeneric("nData", function(x, ...) {
     standardGeneric("nData")
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


