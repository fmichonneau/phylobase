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
  if (length(tip)==0) phy else as(ape::drop.tip(as(phy,"phylo"),tip,...),class(phy))
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
            oldnodelabels <- phy@node.label
            ## need unique labels to match data correctly
            tags <- paste("N",1:nNodes(phy),sep="")
            phy@node.label <- tags
            rownames(phy@node.data) <- phy@node.label
            phytr <- DropTip(phy,tip,trim.internal, subtree, root.edge)
            ## this DROPS data
            phytr@tip.data=phy@tip.data[phytr@tip.label,,drop=FALSE]
            m1  = match(phytr@node.label,tags)
            phytr@node.data=phy@node.data[m1,,drop=FALSE]
            phytr@node.label=oldnodelabels[m1]
            rownames(phytr@node.data)=phytr@node.label
            phytr@node.label=oldnodelabels
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

