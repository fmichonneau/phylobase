setMethod("reorder", signature(x = "phylo"), function(x, order = 'cladewise') {
    x <- ape::reorder.phylo(x,  order)
    x
})
