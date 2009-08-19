## drop.tip4.R

## modified from drop.tip.R (2006-10-12)
##   Remove Tips in a Phylogenetic Tree
## Copyright 2003-2006 Emmanuel Paradis
## This file is part of the R-package `ape'.

## See the file ../COPYING for licensing issues.

setGeneric("prune",function(phy, ...) {
  standardGeneric("prune")
})

## setGeneric("drop.tip") ## if ape has already been loaded


DropTip <- function(phy,tip,...) {
  if (length(tip)==0) {
      phy
  } else if (is(phy,"phylo4d")) {
      ## use extract.tree instead of as() to avoid warning
      as(ape::drop.tip(as(extractTree(phy),"phylo"),tip,...),"phylo4d")
  } else as(ape::drop.tip(as(phy,"phylo"),tip,...),class(phy))
}


## return characters, sorted in NUMERIC order
.chnumsort <- function(x) {
  as.character(sort(as.numeric(x)))
}

setMethod("prune","phylo4",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   ...) {
            DropTip(phy,tip,trim.internal, subtree)
          })

## trace("prune", browser, signature = "phylo4d")
## untrace("prune", signature = "phylo4d")
setMethod("prune", "phylo4d", function(phy, tip, trim.internal=TRUE,
                                       subtree=FALSE, ...) {
    tree <- extractTree(phy)
    phytr <- DropTip(tree, tip, trim.internal, subtree)

    ## create temporary phylo4 object with unique labels
    tmpLbl <- .genlab("n", nTips(phy)+nNodes(phy))
    tmpPhy <- tree
    labels(tmpPhy, "all") <- tmpLbl
    tmpPhytr <- DropTip(tmpPhy, getNode(phy, tip), trim.internal, subtree)

    ## get node numbers to keep
    oldLbl <- labels(tmpPhy, "all")
    newLbl <- labels(tmpPhytr, "all")
    toKeep <- as.numeric(names(oldLbl[oldLbl %in% newLbl]))
    tipToKeep <- toKeep[toKeep %in% nodeId(phy, "tip")]
    nodToKeep <- toKeep[toKeep %in% nodeId(phy, "internal")]

    if(!all(dim(phy@tip.data) == 0)) {
        tipDt <- phy@tip.data[match(tipToKeep, rownames(phy@tip.data)) ,, drop=FALSE]
        tipDt <- tipDt[.chnumsort(rownames(tipDt)) ,, drop=FALSE]
        rownames(tipDt) <- 1:nTips(phytr)
    }
    else
        tipDt <- data.frame(NULL)

    if(!all(dim(phy@node.data) == 0)) {
        nodDt <- phy@node.data[match(nodToKeep, rownames(phy@node.data)) ,, drop=FALSE]
        nodDt <- nodDt[.chnumsort(rownames(nodDt)) ,, drop=FALSE]
        rownames(nodDt) <- 1:nNodes(phytr)
    }
    else
        nodDt <- data.frame(NULL)

    phytr <- phylo4d(phytr, tip.data=tipDt, node.data=nodDt, match.data=FALSE)

    phytr
})

setMethod("prune", "phylo",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   ...) {
            DropTip(phy, tip, trim.internal, subtree)
          })

## setMethod("prune","ANY",
##           function(phy, tip, trim.internal = TRUE, subtree = FALSE,
##                    ,...) {
##             if (class(phy)=="phylo") {
##               ape::prune(phy, tip, trim.internal, subtree)
##               } else stop("no prune method available for",
##                     deparse(substitute(phy)),
##                     "(class",class(phy),")")
##           })

