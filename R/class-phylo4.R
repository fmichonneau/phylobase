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
## Labels constructor
#####################

.createLabels <- function(value, ntips, nnodes, use.names = TRUE,
                          type = c("tip", "internal", "allnode")) {

    type <- match.arg(type)

    ## set up final length of object to return
    lgthRes <- switch(type, tip=ntips, internal=nnodes, allnode=ntips+nnodes)

    ## create NA character vector of node labels
    res <- character(lgthRes)
    is.na(res) <- TRUE

    ## create internal names
    names(res) <- switch(type,
                         tip = 1:ntips,
                         internal = seq(from=ntips+1, length=lgthRes),
                         allnode = 1:(ntips+nnodes))


    ## if no values are provided
    if(missing(value) || is.null(value) || all(is.na(value))) {
        ## tip labels can't be NULL
        if(!identical(type, "internal")) {
            tipLbl <- .genlab("T", ntips)
            res[1:ntips] <- tipLbl
        }
    }

    ## if labels are provided
    else {
        ## check that lengths match
        if(length(value) != lgthRes)
            stop("Number of labels does not match number of nodes.")

        ## check if vector 'value' has name, and if so match with node.label names
        if(use.names && !is.null(names(value))) {
            if(!all(names(value) %in% names(res)))
                stop("Names provided don't match internal labels names.")
            res[match(names(value), names(res))] <- value
        }
        else
            res[1:lgthRes] <- value
    }

    res
}


.createEdge <- function(value, edgeMat, type=c("lengths", "labels"), use.names=TRUE) {
    type <- match.arg(type)

    lgthRes <- nrow(edgeMat)
    res <- switch(type, lengths=numeric(lgthRes), labels=character(lgthRes))
    is.na(res) <- TRUE
    names(res) <- paste(edgeMat[,1], edgeMat[,2], sep="-")

    if(!(missing(value) || is.null(value) || all(is.na(value)))) {
        if(use.names && !is.null(names(value))) {
            if(!all(names(value) %in% names(res)))
                stop("Names provided don't match internal edge labels names.")
            res[match(names(value), names(res))] <- value
        }
        else
            res[1:lgthRes] <- value
    }

    res
}

#####################
## phylo4 constructor
#####################

## generic
setGeneric("phylo4", function(x, ...) { standardGeneric("phylo4")} )

# ape orderings should be allowed for so we can import trees from ape e.g. during subsetting
phylo4_orderings <- c("unknown", "preorder", "postorder", "pruningwise", "cladewise")

## first arg is a matrix
setMethod("phylo4", "matrix",
    function(x, edge.length = NULL, tip.label = NULL, node.label = NULL,
             edge.label = NULL, order="unknown", ...) {

    ## edge
    edge <- x
    mode(edge) <- "integer"
    #if(any(is.na(edge))) stop("NA are not allowed in edge matrix")
    if(ncol(edge) > 2)
        warning("The edge matrix has more than two columns, ",
                "only the first two columns are considered.")
    edge <- as.matrix(edge[, 1:2])
    colnames(edge) <- c("ancestor", "descendant")

    ## number of tips and number of nodes
    ntips <- sum(tabulate(na.omit(edge[, 1])) == 0)
    nnodes <- length(unique(na.omit(c(edge)))) - ntips

    ## edge.length
    edge.length <- .createEdge(value=edge.length, edgeMat=edge, type="lengths", use.names=FALSE)

    ## edge.label
    edge.label <- .createEdge(value=edge.label, edgeMat=edge, type="labels", use.names=FALSE)

    ## tip.label
    tip.label <- .createLabels(value=tip.label, ntips=ntips, nnodes=nnodes,
                               type="tip")

    ## node.label
    node.label <- .createLabels(node.label, ntips=ntips, nnodes=nnodes,
                                type="internal")

    ## fill in the result
    res <- new("phylo4")
    res@edge <- edge
    res@edge.length <- edge.length
    res@Nnode <- nnodes
    res@tip.label <- tip.label
    res@node.label <- node.label
    res@edge.label <- edge.label
    res@order <- order

    ## checkPhylo4 will return a character string if object is
    ##  bad, otherwise TRUE
    if (is.character(checkval <- checkPhylo4(res))) stop(checkval)
    return(res)
})

## first arg is a phylo
setMethod("phylo4", c("phylo"), function(x, check.node.labels=c("keep",
  "drop")){

  check.node.labels <- match.arg(check.node.labels)
  if (check.node.labels == "drop") x$node.label <- NULL
  res <- as(x, "phylo4")

  return(res)
})
