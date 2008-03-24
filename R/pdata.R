## define class for traits
ptypes <- c("multitype","binary","continuous","DNA","RNA","aacid",
            "other","unknown")

setClass("pdata", representation(data="data.frame",
                                 type="factor",
                                 comment="character",
                                 metadata="list"),
         prototype=list(data=data.frame(),type=factor(),
           comment=character(0),metadata=list()))

## pdata constructor
pdata <- function(data,type,comment,metadata) {
  nvar <- ncol(data)
  if (missing(type)) {
    type <- factor(rep("unknown",nvar),levels=ptypes)
  }
  if (length(type)==1) type <- rep(type,length.out=nvar)
  type <- factor(as.character(type),levels=ptypes)
  if (length(comment)==1) comment <- rep(comment,length.out=nvar)
  obj <- new("pdata",data=data,type=type,comment=comment,metadata)
  check_pdata(obj)
  obj
}

check_pdata <- function(object) {
    nvar <- ncol(object@data)
    badlevels <- levels(object@type)[!levels(object@type) %in% ptypes]
    if (length(badlevels)>0)
      stop(paste("bad levels in types:",paste(badlevels,collapse=",")))
    if (length(object@comment)>1 && length(object@comment)!=nvar) {
        stop("wrong number of comments")
    }
    if (length(object@type)>1 && length(object@type)!=nvar) {
        stop("wrong number of types")
    }
}

setMethod("[","pdata",function(x,i, j,...,drop=FALSE) {
  xd <- x@data[i,j,...,drop=drop]
  xd2 <- as.data.frame(xd)
  xd2
})

setMethod("[<-","pdata",function(x,i, j,...,drop=FALSE,value) {
  "[<-"(x@data,i,j,...,drop=drop,value)
})

setGeneric("[[")
setMethod("[[","pdata",
          function(x,i,j,...,exact=NA) {
            x@data[[i,j,...,exact=exact]]
          })

setGeneric("[[<-")
setMethod("[[<-","pdata",
          function(x,i,j,...,exact=NA,value) {
            "[[<-"(x@data,i,j,...,exact=exact,value)
          })

setMethod("$","pdata",function(x,name) {
  x@data[[name]]
})


setMethod("$<-","pdata",function(x,name,value) {
  x@data[[name]] <- value
  x
})

setMethod("plot",signature(x="pdata",y="missing"), function(x,...){
    return(plot(x@data, ...))
}) # end plot phylo4


## od = data.frame(a=1:3,b=4:6)
## z = new("pdata",
##   data=od,type=factor("a","b"),
##   comment=c("",""),metadata=list())

## z[2,]
## z[,"a"]
## z[[2]]

## test conflict resolution error