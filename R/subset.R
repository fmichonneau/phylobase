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
        is.valid.tip <- nodes %in% all.tips
        if (sum(is.valid.tip)<2) {
            stop("mrca must include at least two valid tips")
        }
        mnode <- MRCA(x, nodes[is.valid.tip])
        kept <- descendants(x, mnode)
        dropped <- setdiff(all.tips, kept)
        unknown <- mrca[!is.valid.tip]
    } else if (!is.null(node.subtree)) {
        node <- getNode(x, node.subtree, missing="OK")
        if (length(node)!=1 || !(node %in% nodeId(x, "internal"))) {
            stop("node.subtree must be a single valid internal node")
        }
        kept <- descendants(x, node)
        dropped <- setdiff(all.tips, kept)
        unknown <- numeric(0)
    } else {
        kept <- x@tip.label
        dropped <- numeric(0)
        unknown <- numeric(0)
    }
    if (length(unknown)>0) {
        warning("unknown tips ignored: ", paste(unknown, 
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
## phylo4
setMethod("[", "phylo4", function(x, i, j, ..., drop=FALSE) {
    if (missing(i)) i <- TRUE
    if (is.logical(i)) {
        i <- nodeId(x, "tip")[i]
    } # else pass 'i' straight through to subset method
    subset(x, tips.include=i)
})

## phylo4d
setMethod("[","phylo4d", 
          function(x, i, j, ..., drop=FALSE) {

              if(missing(i)) i <- TRUE
              if(missing(j)) j <- TRUE

              #### data handling
              ## for now handle only tip data - assumes tip names are good row.names
              tab <- tdata(x, type="tip")[i, j, ..., drop=FALSE]
              
              #### tree handling
              tip.include <- match(row.names(tab), x@tip.label)
              tre <- subset(as(x,"phylo4"), tips.include=tip.include)

              ## result
              res <- phylo4d(x=tre, tip.data=tab)
              
              return(res)
          })

## extract the phylo4 part of phylo4d; relies on implicit coerce method
extractTree <- function(from) {
    as(from, "phylo4")
}
