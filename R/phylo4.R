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

setGeneric("tdata<-", function(object,...,value) {
    standardGeneric("tdata<-")
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


