setMethod("print", "phylo4d", printphylo4)

setMethod("show", "phylo4d", function(object) printphylo4(object))

setMethod("tdata", "phylo4d", function(x, which = c("tip", 
    "node", "allnode"), ...) {
    which <- match.arg(which)
    if (which == "allnode") {
        namesmatch <- all(colnames(x@tip.data) == colnames(x@node.data))
        classmatch <- all(sapply(x@tip.data, class) == sapply(x@node.data, 
            class))
        if (!(classmatch && namesmatch)) 
            stop("Node and tip columns do not match, access tip and node data separately")
    }
    switch(which, tip = x@tip.data, node = x@node.data, allnode = rbind(x@tip.data, 
        x@node.data))
})

setMethod("tdata<-", "phylo4d", function(object, which = c("tip", 
    "node", "allnode"), ..., value) {
    which <- match.arg(which)
    if (which == "allnode") {
        namesmatch <- all(colnames(object@tip.data) == colnames(object@node.data))
        classmatch <- all(sapply(object@tip.data, class) == sapply(object@node.data, 
            class))
        if (!(classmatch && namesmatch)) 
            stop("Node and tip columns do not match;",
                 "you should access tip and node data separately")
    }
    switch(which,
           ## FIXME: add checks for matching row names etc ...
           tip = object@tip.data <- value,
           node = object@node.data <- value,
           allnode = stop("for now, must set tip and node data separately"))
    object
})



## Alternative phylo4d summary method, using phylo4 summary
## Marguerite Butler & Peter Cowan
setMethod("summary", "phylo4d", function(object) {
    x <- object
    summary(as(object, "phylo4"))
    tips <- tdata(object, "tip")
    nodes <- tdata(object, "node")
    cat("\nComparative data:\n")
    if (nrow(tips) > 0) {
        cat("\nTips: data.frame with", nTips(object), "taxa and", 
            ncol(tips), "variables \n\n")
        print(summary(tips))
    }
    else {
        cat("\nObject contains no tip data.")
    }
    if (nrow(nodes) > 0) {
        cat("\nNodes: data.frame with", nNodes(object), "internal nodes and", 
            ncol(nodes), "variables \n\n")
        print(summary(nodes))
    }
    else {
        cat("\nObject contains no node data.\n")
    }
})

setMethod("hasNodeData", "phylo4d", function(x) {
    nrow(x@node.data) > 0
})

setMethod("nodeLabels<-", "phylo4d", function(object, ..., 
    value) {
    object@node.label <- value
    rownames(object@node.data) <- value
    object
})

setMethod("labels<-", "phylo4d", function(object, ..., value) {
    object@tip.label <- value
    rownames(object@tip.data) <- value
    object
})

## FIXME: doesn't deal with missing node data
##   (don't even know how that should be done in this case)
setMethod("na.omit", "phylo4d", function(object, ...) {
    tipdata <- tdata(object, "tip")
    na.names <- rownames(tipdata)[!complete.cases(tipdata)]
    prune(object, tip = na.names)
})

setMethod("names", signature(x = "phylo4d"), function(x) {
    temp <- rev(names(attributes(x)))[-1]
    return(rev(temp))
})
