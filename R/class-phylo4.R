setClass("phylo4", 
         representation(edge = "matrix", 
                        edge.length = "numeric", 
                        Nnode = "integer", 
                        node.label = "character", 
                        tip.label = "character", 
                        edge.label = "character", 
                        root.edge = "numeric"), 
         prototype = list(
                        edge = matrix(nrow = 0, ncol = 2, 
                            dimname = list(NULL, c("ancestor", "descendant"))), 
                        edge.length = numeric(0), 
                        Nnode = as.integer(0), 
                        tip.label = character(0), 
                        node.label = character(0), 
                        edge.label = character(0), 
                        root.edge = as.numeric(NA)
                       ), 
         validity = check_phylo4)

#####################
## phylo4 constructor
#####################

phylo4 <- function(edge, edge.length = NULL, tip.label = NULL, node.label = NULL, 
                   edge.label = NULL, root.edge = NULL, ...){
    ## edge
    mode(edge) <- "integer"
    if(any(is.na(edge))) stop("NA are not allowed in edge matrix")
    if(ncol(edge) > 2) warning("the edge matrix has more than two columns")
    edge <- as.matrix(edge[, 1:2])
    colnames(edge) <- c("ancestor", "descendant")

    ## edge.length
    if(!is.null(edge.length)) {
        if(!is.numeric(edge.length)) stop("edge.length is not numeric")
        edge.length <- edge.length
    } else {
        edge.length <- as.numeric(NULL)
    }

    ## tip.label
    ntips <- sum(tabulate(edge[, 1]) == 0)
    if(is.null(tip.label)) {
        tip.label <- .genlab("T", ntips)
    } else {
        if(length(tip.label) != ntips) stop("the tip labels are not consistent with the number of tips")
        tip.label <- as.character(tip.label)
    } 

    ## node.label
    nnodes <- sum(tabulate(edge[, 1]) > 0)
    if(is.null(node.label)) {
        node.label <- .genlab("N", nnodes)
    } else {
        if(length(node.label) != nnodes) stop("the node labels are not consistent with the number of nodes")
    } 

    ## edge.label
    ## an edge is named by the descendant
    if(is.null(edge.label)) {
        edge.label <- paste("E", edge[, 2], sep = "")
    } else {
        if(length(edge.label) != nrow(edge)) stop("the edge labels are not consistent with the number of edges")
    }

    ## root.edge - if no root edge lenth provided, set to a numeric NA
    if(is.null(root.edge)) root.edge <- as.numeric(NA)
    ##if(!is.null(root.edge)) {
    ##    if(!round(root.edge)==root.edge) stop("root.edge must be an integer")
    ##    root.edge <- as.integer(root.edge)
    ##    if(root.edge > nrow(edge)) stop("indicated root.edge do not exist")
    ##} else {
    ##    root.edge <- as.integer(NA)
    ##}

    ## fill in the result
    res <- new("phylo4")
    res@edge <- edge
    res@edge.length <- edge.length
    res@Nnode <- nnodes
    res@tip.label <- tip.label
    res@node.label <- node.label
    res@edge.label <- edge.label
    res@root.edge <- root.edge

    ## check_phylo4 will return a character string if object is
    ##  bad, otherwise TRUE
    if (is.character(checkval <- check_phylo4(res))) stop(checkval)
    return(res)
}

