#' The phylo4 class
#' 
#' Classes for phylogenetic trees
#' 
#' 
#' @name phylo4-class
#' @aliases phylo4_orderings phylo-class phylo4-class
#' @docType class
#' @section Objects from the Class: Phylogenetic tree objects can be created by
#' calls to the \code{\link{phylo4}} constructor function.  Translation
#' functions from other phylogenetic packages are also available. See
#' \code{\link{coerce-methods}}.
#' @author Ben Bolker, Thibaut Jombart
#' @seealso The \code{\link{phylo4}} constructor, the \code{\link{checkPhylo4}}
#' function to check the validity of \code{phylo4} objects. See also the
#' \code{\link{phylo4d}} constructor and the \linkS4class{phylo4d} class.
#' @keywords classes
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


#' Create a phylogenetic tree
#' 
#' \code{phylo4} is a generic constructor that creates a phylogenetic tree
#' object for use in phylobase methods. Phylobase contains functions for input
#' of phylogenetic trees and data, manipulation of these objects including
#' pruning and subsetting, and plotting. The phylobase package also contains
#' translation functions to forms used in other comparative phylogenetic method
#' packages.
#' 
#' The minimum information necessary to create a phylobase tree object is a
#' valid edge matrix. The edge matrix describes the topology of the phylogeny.
#' Each row describes a branch of the phylogeny, with the (descendant) node
#' number in column 2 and its ancestor's node number in column 1. These numbers
#' are used internally and must be unique for each node.
#' 
#' The labels designate either nodes or edges. The vector \code{node.label}
#' names internal nodes, and together with \code{tip.label}, name all nodes in
#' the tree. The vector \code{edge.label} names all branches in the tree. All
#' label vectors are optional, and if they are not given, internally-generated
#' labels will be assigned. The labels, whether user-specified or internally
#' generated, must be unique as they are used to join species data with
#' phylogenetic trees.
#' 
#' @name phylo4-methods
#' @aliases phylo4 phylo4-methods phylo4,matrix-method phylo4,phylo-method
#' @docType methods
#' @param x a matrix of edges or an object of class \code{phylo} (see above)
#' @param edge A numeric, two-column matrix with as many rows as branches in
#' the phylogeny.
#' @param edge.length Edge (branch) length. (Optional)
#' @param tip.label A character vector of species names (names of "tip" nodes).
#' (Optional)
#' @param node.label A character vector of internal node names. (Optional)
#' @param edge.label A character vector of edge (branch) names. (Optional)
#' @param order character: tree ordering (allowable values are listed in
#' \code{phylo4_orderings}, currently "unknown", "preorder" (="cladewise" in
#' \code{ape}), and "postorder", with "cladewise" and "pruningwise" also
#' allowed for compatibility with \code{ape})
#' @param check.node.labels if \code{x} is of class \code{phylo}, either "keep"
#' (the default) or "drop" node labels. This argument is useful if the
#' \code{phylo} object has non-unique node labels.
#' @param annote any additional annotation data to be passed to the new object
#' @note Translation functions are available from many valid tree formats. See
#' \link{coerce-methods}.
#' @section Methods: \describe{ \item{x = "matrix"}{creates a phylobase tree
#' from a matrix of edges}
#' 
#' \item{x = "phylo"}{creates a phylobase tree from an object of class
#' \code{phylo}} }
#' @author phylobase team
#' @seealso \code{\link{coerce-methods}} for translation functions. The
#' \linkS4class{phylo4} class, the \code{\link{formatData}} function to check
#' the validity of \code{phylo4} objects. See also the \code{\link{phylo4d}}
#' constructor, and \linkS4class{phylo4d} class.
#' @keywords classes
#' @examples
#' 
#' # a three species tree:
#' mytree <- phylo4(x=matrix(data=c(4,1, 4,5, 5,2, 5,3, 0,4), ncol=2,
#' byrow=TRUE), tip.label=c("speciesA", "speciesB", "speciesC")) 
#' mytree
#' plot(mytree)
#' 
#' # another way to specify the same tree:
#' mytree <- phylo4(x=cbind(c(4, 4, 5, 5, 0), c(1, 5, 2, 3, 4)),
#' tip.label=c("speciesA", "speciesB", "speciesC")) 
#' 
#' # another way:
#' mytree <- phylo4(x=rbind(c(4, 1), c(4, 5), c(5, 2), c(5, 3), c(0, 4)),
#' tip.label=c("speciesA", "speciesB", "speciesC")) 
#' 
#' # with branch lengths:
#' mytree <- phylo4(x=rbind(c(4, 1), c(4, 5), c(5, 2), c(5, 3), c(0, 4)),
#' tip.label=c("speciesA", "speciesB", "speciesC"), edge.length=c(1, .2,
#' .8, .8, NA))
#' plot(mytree)
#' 
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
    #if(any(is.na(edge))) stop("NA are not allowed in edge matrix") ## taken care by checkTree
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

    ## edge.length (drop elements if all are NA but keep the vector named)
    edge.length <- .createEdge(value=edge.length, edgeMat=edge, type="lengths", use.names=FALSE)
    if (all(is.na(edge.length))) {
        edge.length <- numeric()
        attributes(edge.length) <- list(names=character(0))
    }

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
