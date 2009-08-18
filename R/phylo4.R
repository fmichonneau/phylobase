setGeneric("nNodes", function(x) {
    standardGeneric("nNodes")
})

setGeneric("nEdges", function(x) {
    standardGeneric("nEdges")
})

setGeneric("edges", function(x, order, ...) {
    standardGeneric("edges")
})

setGeneric("rootEdge", function(x, order, ...) {
    standardGeneric("rootEdge")
})

setGeneric("nodeType", function(phy) {
    standardGeneric("nodeType")
})

setGeneric("isRooted", function(x) {
    standardGeneric("isRooted")
})

setGeneric("rootNode", function(x) {
    standardGeneric("rootNode")
})

setGeneric("rootNode<-", function(x, value) {
    standardGeneric("rootNode<-")
})

setGeneric("hasEdgeLength", function(x) {
    standardGeneric("hasEdgeLength")
})

setGeneric("edgeLength", function(x, ...) {
    standardGeneric("edgeLength")
})

setGeneric("edgeLength<-", function(x, ..., value) {
    standardGeneric("edgeLength<-")
})

setGeneric("sumEdgeLength", function(phy, node) {
    standardGeneric("sumEdgeLength")
})

setGeneric("hasNodeLabels", function(x) {
    standardGeneric("hasNodeLabels")
})

setGeneric("hasEdgeLabels", function(x) {
    standardGeneric("hasEdgeLabels")
})

setGeneric("edgeOrder", function(x, ...) {
  standardGeneric("edgeOrder")
})

setGeneric("labels")

setGeneric("labels<-",
           function(object, type, ..., value) {
               standardGeneric("labels<-")
           })

setGeneric("nodeLabels", function(object) {
    standardGeneric("nodeLabels")
})

setGeneric("nodeLabels<-",
           function(object, ..., value) {
               standardGeneric("nodeLabels<-")
           })

setGeneric("tipLabels", function(object) {
    standardGeneric("tipLabels")
})

setGeneric("tipLabels<-",
   function(object, ..., value) {
       standardGeneric("tipLabels<-")
   })

setGeneric("nodeId", function(x, type=c("internal", "tip", "all")) {
    standardGeneric("nodeId")
})

setGeneric("edgeLabels", function(x) {
    standardGeneric("edgeLabels")
})

setGeneric("edgeLabels<-",
           function(object, ..., value) {
               standardGeneric("edgeLabels<-")
           })

setGeneric("print")

setGeneric("head")

setGeneric("tail")

setGeneric("tdata", function(x, ...) {
    standardGeneric("tdata")
})

setGeneric("tdata<-", function(object, ..., value) {
    standardGeneric("tdata<-")
})

setGeneric("addData", function(x, ...) {
    standardGeneric("addData")
})

setGeneric("hasNodeData", function(x) {
    standardGeneric("hasNodeData")
})

setGeneric("na.omit")

setGeneric("reorder")

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


