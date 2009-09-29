################
## subset phylo4
################

setGeneric("subset")
setMethod("subset", "phylo4", function(x, tips.include=NULL,
    tips.exclude=NULL, mrca=NULL, node.subtree=NULL, ...) {
    ##  FIXME: could eliminate NULL and make the test
    ## if (!missing) rather than if (!is.null)
    ## (might have to change next line?)
    if (sum(!sapply(list(tips.include, tips.exclude, mrca,
        node.subtree), is.null))>1) {
        stop("must specify at most one criterion for subsetting")
    }
    #arglist <- list(...)
    #if (length(arglist)>0) {
    #  warning("unused arguments: ",
    #          paste(names(arglist),collapse=","))
    #}
    all.tips <- nodeId(x, "tip")
    if (!is.null(tips.include)) {
        nodes <- getNode(x, tips.include, missing="OK")
        is.valid.tip <- nodes %in% all.tips
        kept <- nodes[is.valid.tip]
        dropped <- setdiff(all.tips, kept)
        unknown <- tips.include[!is.valid.tip]
    } else if (!is.null(tips.exclude)) {
        nodes <- getNode(x, tips.exclude, missing="OK")
        is.valid.tip <- nodes %in% all.tips
        dropped <- nodes[is.valid.tip]
        kept <- setdiff(all.tips, dropped)
        unknown <- tips.exclude[!is.valid.tip]
    } else if (!is.null(mrca)) {
        nodes <- getNode(x, mrca, missing="OK")
        is.valid.node <- nodes %in% nodeId(x, "all")
        mnode <- MRCA(x, nodes[is.valid.node])
        if (length(mnode)!=1) {
            stop("mrca must include at least one valid node")
        }
        kept <- descendants(x, mnode)
        dropped <- setdiff(all.tips, kept)
        unknown <- mrca[!is.valid.node]
    } else if (!is.null(node.subtree)) {
        node <- getNode(x, node.subtree, missing="OK")
        if (length(node)!=1 || !(node %in% nodeId(x, "internal"))) {
            stop("node.subtree must be a single valid internal node")
        }
        kept <- descendants(x, node)
        dropped <- setdiff(all.tips, kept)
        unknown <- numeric(0)
    } else {
        kept <- getNode(x, nodeId(x, "tip"))
        dropped <- numeric(0)
        unknown <- numeric(0)
    }
    if (length(unknown)>0) {
        warning("invalid nodes ignored: ", paste(unknown, 
            collapse=", "))
    }
    if (length(kept)<2) {
        stop("0 or 1 tips would remain after subsetting")
    }
    if (length(dropped)==0) return(x)
    return(prune(x, dropped, ...))
})

###############
# '[' operator
###############

## Consider using some combination of these for stricter argument
## checking? Not implementing now because extra arguments are just
## ignored, which is fairly common S4 method behavior:
## * in "[" methods for phylo4:
##    if (nargs()>2) stop("unused arguments")
## * in "[" methods for both phylo4 and phylo4d:
##    if (!missing(...)) stop("unused argument(s)")

## phylo4 '[' methods
setMethod("[", signature(x="phylo4", i="character", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=i)
})
setMethod("[", signature(x="phylo4", i="numeric", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=i)
})
setMethod("[", signature(x="phylo4", i="logical", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=nodeId(x, "tip")[i])
})
setMethod("[", signature(x="phylo4", i="missing", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    x
})
## phylo4d '[' methods
setMethod("[", signature(x="phylo4d", i="ANY", j="character",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})
setMethod("[", signature(x="phylo4d", i="ANY", j="numeric",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})
setMethod("[", signature(x="phylo4d", i="ANY", j="logical",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})
## borrow from Matrix package approach of trapping invalid usage
setMethod("[", signature(x="phylo4", i="ANY", j="ANY", drop="ANY"),
    function(x, i, j, ..., drop) {
    stop("invalid argument(s)")
})


## extract the phylo4 part of phylo4d; relies on implicit coerce method
extractTree <- function(from) {
    as(from, "phylo4")
}
