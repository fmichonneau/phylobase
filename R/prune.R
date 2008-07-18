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
      as(ape::drop.tip(as(extract.tree(phy),"phylo"),tip,...),"phylo4d")
  } else as(ape::drop.tip(as(phy,"phylo"),tip,...),class(phy))
}


setMethod("prune","phylo4",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0,...) {
            DropTip(phy,tip,trim.internal, subtree, root.edge)
          })

## trace("prune", browser, signature = "phylo4d")
## untrace("prune", signature = "phylo4d")
setMethod("prune","phylo4d",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0,...) {
            ## need unique labels to match data correctly
            oldnodelabels <- phy@node.label
            nodetags <- .genlab("N",nNodes(phy))
            phy@node.label <- nodetags
            oldtiplabels <- phy@tip.label
            phytr <- DropTip(phy,tip,trim.internal, subtree, root.edge)
            ## this DROPS data
            ntr = match(phytr@node.label,nodetags)
            ttr = match(phytr@tip.label,oldtiplabels)
            phytr@node.label <- oldnodelabels[ntr]
            phytr@tip.label <- oldtiplabels[ttr]
            phytr@node.data <- phy@node.data[ntr,,drop=FALSE]
            phytr@tip.data <- phy@tip.data[ttr,,drop=FALSE]            
            phytr
          })

setMethod("prune","phylo",
          function(phy, tip, trim.internal = TRUE, subtree = FALSE,
                   root.edge = 0,...) {
            DropTip(phy,tip,trim.internal, subtree, root.edge)
          })

## setMethod("prune","ANY",
##           function(phy, tip, trim.internal = TRUE, subtree = FALSE,
##                    root.edge = 0,...) {
##             if (class(phy)=="phylo") {
##               ape::prune(phy, tip, trim.internal, subtree, root.edge)
##               } else stop("no prune method available for",
##                     deparse(substitute(phy)),
##                     "(class",class(phy),")")
##           })

