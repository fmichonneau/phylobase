################
## subset phylo4
################
##

## 

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
                  return(prune(x,x@tip.label[!(x@tip.label %in% allDescend(x,node.subtree))]))
              }
              
              if (!is.null(mrca)) {
                  mnode <- MRCA(x,mrca)
                  return(prune(x,x@tip.label[!(x@tip.label %in% allDescend(x,mnode))]))
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


