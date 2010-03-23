setClass("phylo4",
         representation(edge = "matrix",
                        edge.length = "numeric",
                        label = "character",
                        edge.label = "character",
                        order = "character",
                        annote = "list"),
         prototype = list(
                        edge = matrix(nrow = 0, ncol = 2,
                            dimname = list(NULL, c("ancestor", "descendant"))),
                        edge.length = numeric(0),
                        label = character(0),
                        edge.label = character(0),
                        order = "unknown",
                        annote = list()
                       ),
         validity = checkPhylo4)

#####################
## Labels constructor
#####################

.createLabels <- function(value, ntips, nnodes, use.names = TRUE,
                          type = c("all", "tip", "internal")) {

    type <- match.arg(type)

    ## set up final length of object to return
    lgthRes <- switch(type, tip=ntips, internal=nnodes, all=ntips+nnodes)

    ## create NA character vector of node labels
    res <- character(lgthRes)
    is.na(res) <- TRUE

    ## create internal names
    names(res) <- switch(type,
                         tip = 1:ntips,
                         internal = seq(from=ntips+1, length=lgthRes),
                         all = 1:(ntips+nnodes))

    ## Convert empty labels to NA
    value[!nzchar(value)] <- NA

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
             edge.label = NULL, order="unknown", annote = list()) {

    ## edge
    edge <- x
    mode(edge) <- "integer"
    #if(any(is.na(edge))) stop("NA are not allowed in edge matrix")
    if(ncol(edge) > 2)
        warning("The edge matrix has more than two columns, ",
                "only the first two columns are considered.")
    edge <- as.matrix(edge[, 1:2])
    colnames(edge) <- c("ancestor", "descendant")

    ## create new phylo4 object and insert edge matrix
    res <- new("phylo4")
    res@edge <- edge

    ## get number of tips and number of nodes
    ## (these accessors work fine now that edge matrix exists)
    ntips <- nTips(res)
    nnodes <- nNodes(res)

    ## edge.length (drop elements if all are NA)
    edge.length <- .createEdge(value=edge.length, edgeMat=edge, type="lengths", use.names=FALSE)
    if (all(is.na(edge.length))) edge.length <- numeric()

    ## edge.label (drop NA elements)
    edge.label <- .createEdge(value=edge.label, edgeMat=edge, type="labels", use.names=FALSE)
    edge.label <- edge.label[!is.na(edge.label)]

    ## tip.label (leave NA elements; let checkTree complain about it)
    tip.label <- .createLabels(value=tip.label, ntips=ntips, nnodes=nnodes,
                               type="tip")

    ## node.label (drop NA elements)
    node.label <- .createLabels(node.label, ntips=ntips, nnodes=nnodes,
                                type="internal")
    node.label <- node.label[!is.na(node.label)]

    ## populate the slots
    res@edge.length <- edge.length
    res@label <- c(tip.label, node.label)
    res@edge.label <- edge.label
    res@order <- order
    res@annote <- annote

    ## checkPhylo4 will return a character string if object is
    ##  bad, otherwise TRUE
    if (is.character(checkval <- checkPhylo4(res))) stop(checkval)
    return(res)
})

## first arg is a phylo
setMethod("phylo4", c("phylo"), function(x, check.node.labels=c("keep",
  "drop"), annote=list()){

  check.node.labels <- match.arg(check.node.labels)
  if (check.node.labels == "drop") x$node.label <- NULL
  res <- as(x, "phylo4")
  #TODO?: make default annote arg NULL, and only assign if !is.null;
  # then update phylo4d methods accordingly (same thing with metadata?)
  res@annote <- annote

  return(res)
})
