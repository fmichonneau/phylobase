

setMethod("isRooted", signature(x="phylo4"),
 function(x) {
    ## hack to avoid failure on an empty object
    if(nTips(x) == 0) return(FALSE)
    any(edges(x)[, 1] == 0)
})

setMethod("rootNode", signature(x="phylo4"),
 function(x) {
    if (!isRooted(x))
        return(NA)
    unname(edges(x)[which(edges(x)[, 1] == 0), 2])
})

setReplaceMethod("rootNode", signature(x="phylo4"),
 function(x, value) {
    stop("Root node replacement not implemented yet")
})

