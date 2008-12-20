#######################################################
## Importing from ape
setAs("phylo", "phylo4", function(from, to) {
    #fixme SWK kludgy fix may not work well with unrooted trees
    if (is.rooted(from)) {
        root.edge <- as.numeric(setdiff(unique(from$edge[,1]), unique(from$edge[,2])))
        from$edge <- rbind(from$edge,c(NA,root.edge))
        if (!is.null(from$edge.length)) {
            if (is.null(from$root.edge)) {
                from$edge.length <- c(from$edge.length,as.numeric(NA))
            }
            else {
                from$edge.length <- c(from$edge.length,from$root.edge)
            }
        }
        if (!is.null(from$edge.label)) {
            from$edge.label <- c(from$edge.label,paste("E",root.edge,sep=""))
        }
    }
    newobj <- phylo4(from$edge, from$edge.length, from$tip.label,
        node.label = from$node.label, edge.label = from$edge.label)
    attribs = attributes(from)
    attribs$names <- NULL
    knownattr <- c("logLik", "order", "origin", "para", "xi")
    known <- names(attribs)[names(attribs) %in% knownattr]
    unknown <- names(attribs)[!names(attribs) %in% c(knownattr, "class", "names")]
    if (length(unknown) > 0) {
        warning(paste("unknown attributes ignored: ", unknown, collapse = " "))
    }
    for (i in known) attr(newobj, i) <- attr(from, i)
    newobj
})

setAs("phylo", "phylo4d", function(from, to) {
    phylo4d(as(from, "phylo4"), tip.data = data.frame())
})

setAs("multiPhylo", "multiPhylo4", function(from, to) {
    y <- lapply(as, from@phylolist, to = "phylo")
    names(y) <- from@tree.names
    if (nrow(from@tip.data) > 0)
        warning("discarded tip data")
    class(y) <- "multiPhylo"
    y
})

#######################################################
## Exporting to ape
setAs("phylo4", "phylo", function(from, to) {
  if (inherits(from,"phylo4d"))
    warning("losing data while coercing phylo4d to phylo")
  brlen <- from@edge.length
  if (isRooted(from)) brlen <- brlen[nodeId(from,"all")!=rootNode(from)]
  edgemat <- na.omit(from@edge)
  y <- list(edge = na.omit(from@edge), edge.length = brlen,
            Nnode = from@Nnode, tip.label = from@tip.label,
            node.label = from@node.label)
  class(y) <- "phylo"
  if (length(y$edge.length) == 0)
    y$edge.length <- NULL
  if (length(y$node.label) == 0)
    y$node.label <- NULL
  if (isRooted(from)) {
    root.edge <- brlen[nodeId(from,"all")==rootNode(from)]
    if (!is.na(root.edge)) y$root.edge <- root.edge
  }
  y
})

## BMB: redundant????  
## setAs("phylo4d", "phylo", function(from, to) {
##     y <- list(edge = from@edge, edge.length = from@edge.length,
##         Nnode = from@Nnode, tip.label = from@tip.label)
##     class(y) <- "phylo"
##     if (length(y$edge.length) == 0)
##         y$edge.length <- NULL
##     if (length(y$node.label) == 0)
##         y$node.label <- NULL
##     #if (!is.na(from@root.edge))
##     #    y$root.edge <- from@root.edge
##    warning("losing data while coercing phylo4d to phylo")
##    y
##})

setAs("multiPhylo4", "multiPhylo", function(from, to) {
    newobj <- new("multiPhylo4", phylolist = lapply(from,
        as, to = "phylo4"))
})

setAs("multiPhylo4d", "multiPhylo", function(from, to) {
    newobj <- new("multiPhylo4d", phylolist = lapply(from,
        as, to = "phylo4"), tree.names = names(from), tip.data = data.frame())
})

#######################################################
## Exporting to ade4
setAs("phylo4", "phylog", function(from, to) {
    if (!require(ade4))
        stop("the ade4 package is required")
    x <- as(from, "phylo")
    xstring <- write.tree(x, file = "")
    newick2phylog(xstring)
})

#######################################################
## Exporting to dataframe
setAs(from = "phylo4", to = "data.frame", def = function(from) {
    if (is.character(checkval <- check_phylo4(from))) # check the phylo4
        stop(checkval)
    x <- from
    if (isRooted(x)) {
        E <- edges(x) # E: matrix of edges
        ancestor <- E[, 1]
        node <- E[, 2]
        root <- which(is.na(ancestor))
        int.node <-  node[(node %in% ancestor)]
        tip <- node[!(node %in% ancestor)]
        n.tip <- length(tip)
        n.int <- length(int.node)
        ## node <- c(root, node) # doesn't fit the ordering: root, other internal nodes, tips
        #node <- c(int.node, tip)
        ## retrieve the ancestor of each node
        #idx <- match(node, E[, 2]) # new ordering of the descendants/edges
        ## if (length(ancestor)>0) ancestor <- c(NA, ancestor)
        #ancestor <- E[idx, 1]
        ## branch.length <- c(x@root.edge, x@edge.length) # root.edge is not an edge length
        branch.length <- edgeLength(x)#[idx]
        if (is.null(edgeLength(x))) {
            branch.length <- rep(NA, length(node))
        }
        ## node and tip labels ##
        ## beware: they cannot be NULL
        ## there are always tip labels (or check_phylo4 complains)
        ## there may not be node labels (character(0))
        label <- labels(x,which="all")[node]
        node.type <- nodeType(x)[node]
        return(data.frame(label, node, ancestor, branch.length,
            node.type,stringsAsFactors=FALSE))
    }
    else {
        E <- edges(x) # E: matrix of edges
        node <- unique(c(E))
        ancestor <- E[, 1][node]
        #orphan <- setdiff(E[,1],E[,2])
        branch.length <- edgeLength(x)[node]
        if (is.null(edgeLength(x))) {
            branch.length <- rep(NA, length(node))
        }
        ## node and tip labels ##
        ## beware: they cannot be NULL
        ## there are always tip labels (or check_phylo4 complains)
        ## there may not be node labels (character(0))
        label <- labels(x,which="all")[node]
        node.type <- nodeType(x)[node]
        return(data.frame(label, node, ancestor, branch.length,
            node.type,stringsAsFactors=FALSE))        
    }
})

setAs(from = "phylo4d", to = "data.frame", function(from) {
    tree <- as(from, "phylo4") # get tree
    t_df <- as(tree, "data.frame") # convert to data.frame
    dat <- tdata(from, "allnode", label.type="column") # get data
    tdat <- cbind(t_df,dat[,-1,drop=FALSE])
    #tdat <- dat[,-1,drop=FALSE]
    return(tdat)
})
