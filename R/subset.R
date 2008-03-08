################
## subset phylo4
################

setGeneric("subset")
setMethod("subset", "phylo4",
          function(x,tips.include=NULL,tips.exclude=NULL,
                   mrca=NULL,node.subtree=NULL,...) {
              ##  FIXME: could do eliminate NULL and make the test
              ## if (!missing) rather than if (!is.null)
              if (!is.null(tips.include)) {
                  if (is.numeric(tips.include)) {
                      tips.include <- x@tip.label[tips.include]
                  }
                  return(prune(x,x@tip.label[!(x@tip.label %in% tips.include)]))
              }
              
              if (!is.null(tips.exclude)) {
                  return(prune(x,tips.exclude))
              }
              
              if (!is.null(node.subtree)) {
                  return(prune(x,x@tip.label[!(x@tip.label %in% names(descendants(x,node.subtree)))]))
              }
              
              if (!is.null(mrca)) {
                  mnode <- MRCA(x,mrca)
                  return(prune(x,x@tip.label[!(x@tip.label %in% names(descendants(x,mnode)))]))
              }
              arglist <- list(...)
              if (length(arglist)>0) {
                  warning("unused arguments: ",
                          paste(names(arglist),collapse=","))
              }
              return(x)
          })


setMethod("subset", "phylo", function(x,...) {
    x <- as(x,"phylo4")
    res <- subset(x,...)
    return(as(res,"phylo"))
})




###############
# '[' operator
###############
## phylo4
setMethod("[","phylo4",
          function(x, i) {

              if(missing(i)) i <- TRUE

              oldlab <- labels(x)
              if(is.character(i)){
                  newlab <- i
              } else {
                  newlab <- oldlab[i]
              }
              tip.include <- match(newlab, oldlab)
              res <- subset(x, tips.include=tip.include)

              return(res)
          })




## phylo4d
setMethod("[","phylo4d", 
          function(x, i, j, ..., drop=FALSE) {

              if(missing(i)) i <- TRUE
              if(missing(j)) j <- TRUE

              #### data handling
              ## for now handle only tip data
              tab <- tdata(x, which="tip")[i, j, ...,drop=FALSE]
              oldtabnames <- row.names(tdata(x,which="tip"))
              
              #### tree handling
              tip.include <- match(row.names(tab), oldtabnames)
              tre <- subset(as(x,"phylo4"), tips.include=tip.include)

              ## result
              res <- phylo4d(x=tre, tip.data=tab)
              
              return(res)
          })


## coerce phylo4d to phylo4 -- on purpose, so no warning

extract.tree <- function(from) {
    phylo4(edge = from@edge, edge.length = from@edge.length, 
        Nnode = from@Nnode, tip.label = from@tip.label)
}
