setClass("phylo4",
         representation(edge = "matrix",
                        edge.length = "numeric",
                        Nnode = "integer",
                        node.label = "character",
                        tip.label = "character",
                        edge.label = "character",
                        order = "character"),
         prototype = list(
                        edge = matrix(nrow = 0, ncol = 2,
                            dimname = list(NULL, c("ancestor", "descendant"))),
                        edge.length = numeric(0),
                        Nnode = as.integer(0),
                        tip.label = character(0),
                        node.label = character(0),
                        edge.label = character(0),
                        order = "unknown"
                       ),
         validity = checkPhylo4)

#####################
## phylo4 constructor
#####################

# ape orderings should be allowed for so we can import trees from ape e.g. during subsetting
phylo4_orderings <- c("unknown", "preorder", "postorder", "pruningwise", "cladewise")

phylo4 <- function(edge, edge.length = NULL, tip.label = NULL, node.label = NULL, edge.label = NULL, order="unknown", ...){

    ## edge
    mode(edge) <- "integer"
    #if(any(is.na(edge))) stop("NA are not allowed in edge matrix")
    if(ncol(edge) > 2) warning("the edge matrix has more than two columns")
    edge <- as.matrix(edge[, 1:2])
    colnames(edge) <- c("ancestor", "descendant")

    ## edge.length
    if(!is.null(edge.length)) {
        if(!is.numeric(edge.length)) stop("edge.length is not numeric")
        edge.length <- edge.length
    } else {
        edge.length <- numeric(0)
    }

    if(length(edge.length) > 0) {
        if(length(edge.length) != nrow(edge))
            stop("The number of edge lengths is different from the number of edges.")
        ## FM - 2009-04-19
        ## edge.length is named according to the nodes the edge links together
        ## (ancestor-descendant). This should allow more robust edge/edge.length
        ## association and limit the problems associated with reordering trees.
        names(edge.length) <- paste(edge[,1], edge[,2], sep="-")
    }

    ## tip.label
    ntips <- sum(tabulate(na.omit(edge[, 1])) == 0)
    if(is.null(tip.label)) {
        tip.label <- .genlab("T", ntips)
    } else {
        if(length(tip.label) != ntips)
            stop("the tip labels are not consistent with the number of tips")
        tip.label <- as.character(tip.label)
    }
    names(tip.label) <- seq(along=tip.label)

    ## node.label for internal nodes
    nnodes <- length(unique(na.omit(c(edge)))) - ntips

    if(is.null(node.label)) {
      node.label <- character(0) ## empty node labels
    }
    else {
        if(length(node.label)>0 && length(node.label) != nnodes)
            stop("number of node labels is not consistent with the number of nodes")
    }
    names(node.label) <- seq(from=ntips+1, along=node.label)


    ## edge.label
    if(is.null(edge.label)) {
      edge.label <- character(0)
    } else if (length(edge.label)>0 && length(edge.label) != nrow(edge))
      stop("number of edge labels is not consistent with the number of edges")


    ## fill in the result
    res <- new("phylo4")
    res@edge <- edge
    res@edge.length <- edge.length
    res@Nnode <- nnodes
    res@tip.label <- tip.label
    res@edge.label <- edge.label
    res@order <- order

    ## Tweak to deal with numerical values returned as node labels (ie. MrBayes)
    if(!all(is.na(node.label)) && any(nchar(node.label) > 0) &&
            !length(grep("[a-zA-Z]", node.label))) {
        warning("All node labels are numeric values and converted as data.")
        res@node.label <- character(0)
        node.data <- node.label

        node.data[!nzchar(node.data)] <- NA

        node.label <- character(nnodes)
        is.na(node.label) <- TRUE

        node.data <- data.frame(labelValues=as.numeric(node.data))
        res@node.label <- node.label

        res <- phylo4d(res, node.data=node.data, use.node.names=FALSE)
        if(is.character(checkval <- checkPhylo4(res))) stop(checkval)

        return(res)

    }
    else {
        res@node.label <- node.label
        ## checkPhylo4 will return a character string if object is
        ##  bad, otherwise TRUE
        if (is.character(checkval <- checkPhylo4(res))) stop(checkval)
    }

    return(res)
}

