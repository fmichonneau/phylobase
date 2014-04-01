## classes for holding multiple tree objects

#' multiPhylo4 and extended classes
#' 
#' Classes for lists of phylogenetic trees.  These classes and methods are
#' planned for a future version of \code{phylobase}.
#' 
#' 
#' @name multiPhylo-class
#' @aliases multiPhylo-class multiPhylo4-class multiPhylo4d-class tbind
#' @docType class
#' @keywords classes
#' @export
#' @include class-multiphylo4.R
setClass("multiPhylo4", representation(phylolist = "list", 
    tree.names = "character"), prototype = list(phylolist = list(), 
    tree.names = character(0)))

setClass("multiPhylo4d", representation(tip.data = "data.frame"), 
    contains = "multiPhylo4")

setMethod("initialize", "multiPhylo4", function(.Object, ...) {
    stop("multiPhylo and multiphylo4d not yet implemented", 
         "Try using a list of phylo4(d) objects and lapply().")
})
