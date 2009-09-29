## classes for holding multiple tree objects

setClass("multiPhylo4", representation(phylolist = "list", 
    tree.names = "character"), prototype = list(phylolist = list(), 
    tree.names = character(0)))

setClass("multiPhylo4d", representation(tip.data = "data.frame"), 
    contains = "multiPhylo4")

setMethod("initialize", "multiPhylo4", function(.Object, ...) {
    stop('multiPhylo and multiphylo4d not yet implemented.  
    Try using a list of phylo4(d) objects and lapply()')
})
