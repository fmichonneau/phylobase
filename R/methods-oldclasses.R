setMethod("reorder", signature(x = "phylo"), function(x, order = 'cladewise') {
    x <- reorder.phylo(x,  order)
    x
})
