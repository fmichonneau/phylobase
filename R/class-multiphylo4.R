## classes for holding multiple tree objects

setClass("multiPhylo4", representation(phylolist = "list", 
    tree.names = "character"), prototype = list(phylolist = list(), 
    tree.names = character(0)))

setClass("multiPhylo4d", representation(tip.data = "data.frame"), 
    contains = "multiPhylo4")
