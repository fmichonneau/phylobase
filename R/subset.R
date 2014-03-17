################
## subset phylo4
################

#' Methods for creating subsets of phylogenies
#' 
#' Methods for creating subsets of phylogenies, based on pruning a tree to
#' include or exclude a set of terminal taxa, to include all descendants of the
#' MRCA of multiple taxa, or to return a subtree rooted at a given node.
#' 
#' The \code{subset} methods must be called using no more than one of the four
#' main subsetting criteria arguments (\code{tips.include},
#' \code{tips.exclude}, \code{mrca}, or \code{node.subtree}).  Each of these
#' arguments can be either character or numeric.  In the first case, they are
#' treated as node labels; in the second case, they are treated as node
#' numbers.  For the first two arguments, any supplied tips not found in the
#' tree (\code{tipLabels(x)}) will be ignored, with a warning.  Similarly, for
#' the \code{mrca} argument, any supplied tips or internal nodes not found in
#' the tree will be ignored, with a warning.  For the \code{node.subtree}
#' argument, failure to provide a single, valid internal node will result in an
#' error.
#' 
#' Although \code{prune} is mainly intended as the workhorse function called by
#' \code{subset}, it may also be called directly.  In general it should be
#' equivalent to the \code{tips.exclude} form of \code{subset} (although
#' perhaps with less up-front error checking).
#' 
#' The "[" operator, when used as \code{x[i]}, is similar to the
#' \code{tips.include} form of \code{subset}.  However, the indices used with
#' this operator can also be logical, in which case the corresponding tips are
#' assumed to be ordered as in \code{nodeId(x, "tip")}, and recycling rules
#' will apply (just like with a vector or a matrix).  With a
#' \linkS4class{phylo4d} object 'x', \code{x[i,j]} creates a subset of \code{x}
#' taking \code{i} for a tip index and \code{j} for the index of data variables
#' in \code{tdata(geospiza, "all")}.  Note that the second index is optional:
#' \code{x[i, TRUE]}, \code{x[i,]}, and \code{x[i]} are all equivalent.
#' 
#' Regardless of which approach to subsetting is used, the argument values must
#' be such that at least two tips are retained.
#' 
#' If the most recent common ancestor of the retained tips is not the original
#' root node, then the root node of the subset tree will be a descendant of the
#' original root.  For rooted trees with non-NA root edge length, this has
#' implications for the new root edge length.  In particular, the new length
#' will be the summed edge length from the new root node back to the original
#' root (including the original root edge).  As an alternative, see the
#' examples for a way to determine the length of the edge that was immediately
#' ancestral to the new root node in the original tree.
#' 
#' Note that the correspondance between nodes and labels (and data in the case
#' of \linkS4class{phylo4d}) will be retained after all forms of subsetting.
#' Beware, however, that the node numbers (IDs) will likely be altered to
#' reflect the new tree topology, and therefore cannot be compared directly
#' between the original tree and the subset tree.
#' 
#' @name subset-methods
#' @aliases subset-methods subset,phylo4-method subset,phylo4d-method prune
#' prune-methods prune,phylo4-method prune,phylo4d-method [-methods
#' [,phylo4-method [,phylo4,character,missing,missing-method
#' [,phylo4,numeric,missing,missing-method
#' [,phylo4,logical,missing,missing-method
#' [,phylo4,missing,missing,missing-method [,phylo4d-method
#' [,phylo4d,ANY,character,missing-method [,phylo4d,ANY,numeric,missing-method
#' [,phylo4d,ANY,logical,missing-method [,phylo4d,ANY,missing,missing-method
#' [,phylo4,ANY,ANY,ANY-method na.omit,phylo4d-method
#' @docType methods
#' @param x an object of class \code{"phylo4"} or \code{"phylo4d"}
#' @param tips.include A vector of tips to include in the subset tree
#' @param tips.exclude A vector of tips to exclude from the subset tree
#' @param mrca A vector of nodes for determining the most recent common
#' ancestor, which is then used as the root of the subset tree
#' @param node.subtree A single internal node specifying the root of the subset
#' tree
#' @param trim.internal A logical specifying whether to remove internal nodes
#' that no longer have tip descendants in the subset tree
#' @param i (\code{[} method) An index vector indicating tips to include
#' @param j (\code{[} method, phylo4d only) An index vector indicating columns
#' of node/tip data to include
#' @param \dots additional arguments to be passed to other methods
#' @return an object of class \code{"phylo4"} or \code{"phylo4d"}
#' @section Methods: \describe{ \item{x = "phylo4"}{subset tree} \item{x =
#' "phylo4d"}{subset tree and corresponding node and tip data} }
#' @author Jim Regetz \email{regetz@@nceas.ucsb.edu}\cr Steven Kembel
#' \email{skembel@@berkeley.edu}\cr Damien de Vienne
#' \email{damien.de-vienne@@u-psud.fr}\cr Thibaut Jombart
#' \email{jombart@@biomserv.univ-lyon1.fr}
#' @keywords manip methods
#' @examples
#' 
#' data(geospiza)
#' nodeLabels(geospiza) <- paste("N", nodeId(geospiza, "internal"), sep="")
#' geotree <- extractTree(geospiza)
#' 
#' ## "subset" examples
#' tips <- c("difficilis", "fortis", "fuliginosa", "fusca", "olivacea",
#'     "pallida", "parvulus", "scandens")
#' plot(subset(geotree, tips.include=tips))
#' plot(subset(geotree, tips.include=tips, trim.internal=FALSE))
#' plot(subset(geotree, tips.exclude="scandens"))
#' plot(subset(geotree, mrca=c("scandens","fortis","pauper")))
#' plot(subset(geotree, node.subtree=18))
#' 
#' ## "prune" examples (equivalent to subset using tips.exclude)
#' plot(prune(geotree, tips))
#' 
#' ## "[" examples (equivalent to subset using tips.include)
#' plot(geotree[c(1:6,14)])
#' plot(geospiza[c(1:6,14)])
#' 
#' ## for phylo4d, subset both tips and data columns
#' geospiza[c(1:6,14), c("wingL", "beakD")]
#' 
#' ## note handling of root edge length:
#' edgeLength(geotree)['0-15'] <- 0.1
#' geotree2 <- geotree[1:2]
#' ## in subset tree, edge of new root extends back to the original root
#' edgeLength(geotree2)['0-3']
#' ## edge length immediately ancestral to this node in the original tree
#' edgeLength(geotree, MRCA(geotree, tipLabels(geotree2)))
setGeneric("subset")
setMethod("subset", "phylo4", function(x, tips.include=NULL,
    tips.exclude=NULL, mrca=NULL, node.subtree=NULL, ...) {
    ##  FIXME: could eliminate NULL and make the test
    ## if (!missing) rather than if (!is.null)
    ## (might have to change next line?)
    if (sum(!sapply(list(tips.include, tips.exclude, mrca,
        node.subtree), is.null))>1) {
        stop("must specify at most one criterion for subsetting")
    }
    #arglist <- list(...)
    #if (length(arglist)>0) {
    #  warning("unused arguments: ",
    #          paste(names(arglist),collapse=","))
    #}
    all.tips <- nodeId(x, "tip")
    if (!is.null(tips.include)) {
        nodes <- getNode(x, tips.include, missing="OK")
        is.valid.tip <- nodes %in% all.tips
        kept <- nodes[is.valid.tip]
        dropped <- setdiff(all.tips, kept)
        unknown <- tips.include[!is.valid.tip]
    } else if (!is.null(tips.exclude)) {
        nodes <- getNode(x, tips.exclude, missing="OK")
        is.valid.tip <- nodes %in% all.tips
        dropped <- nodes[is.valid.tip]
        kept <- setdiff(all.tips, dropped)
        unknown <- tips.exclude[!is.valid.tip]
    } else if (!is.null(mrca)) {
        nodes <- getNode(x, mrca, missing="OK")
        is.valid.node <- nodes %in% nodeId(x, "all")
        mnode <- MRCA(x, nodes[is.valid.node])
        if (length(mnode)!=1) {
            stop("mrca must include at least one valid node")
        }
        kept <- descendants(x, mnode)
        dropped <- setdiff(all.tips, kept)
        unknown <- mrca[!is.valid.node]
    } else if (!is.null(node.subtree)) {
        node <- getNode(x, node.subtree, missing="OK")
        if (length(node)!=1 || !(node %in% nodeId(x, "internal"))) {
            stop("node.subtree must be a single valid internal node")
        }
        kept <- descendants(x, node)
        dropped <- setdiff(all.tips, kept)
        unknown <- numeric(0)
    } else {
        kept <- getNode(x, nodeId(x, "tip"))
        dropped <- numeric(0)
        unknown <- numeric(0)
    }
    if (length(unknown)>0) {
        warning("invalid nodes ignored: ", paste(unknown, 
            collapse=", "))
    }
    if (length(kept)<2) {
        stop("0 or 1 tips would remain after subsetting")
    }
    if (length(dropped)==0) return(x)
    return(prune(x, dropped, ...))
})

###############
# '[' operator
###############

## Consider using some combination of these for stricter argument
## checking? Not implementing now because extra arguments are just
## ignored, which is fairly common S4 method behavior:
## * in "[" methods for phylo4:
##    if (nargs()>2) stop("unused arguments")
## * in "[" methods for both phylo4 and phylo4d:
##    if (!missing(...)) stop("unused argument(s)")

## phylo4 '[' methods
setMethod("[", signature(x="phylo4", i="character", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=i)
})
setMethod("[", signature(x="phylo4", i="numeric", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=i)
})
setMethod("[", signature(x="phylo4", i="logical", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    subset(x, tips.include=nodeId(x, "tip")[i])
})
setMethod("[", signature(x="phylo4", i="missing", j="missing",
    drop="missing"), function(x, i, j, ..., drop) {
    x
})
## phylo4d '[' methods
setMethod("[", signature(x="phylo4d", i="ANY", j="character",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})
setMethod("[", signature(x="phylo4d", i="ANY", j="numeric",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})
setMethod("[", signature(x="phylo4d", i="ANY", j="logical",
    drop="missing"), function(x, i, j, ..., drop) {
    if (!missing(i)) x <- x[i]
    tdata(x, type="all") <- tdata(x, type="all")[j]
    return(x)
})
## borrow from Matrix package approach of trapping invalid usage
setMethod("[", signature(x="phylo4", i="ANY", j="ANY", drop="ANY"),
    function(x, i, j, ..., drop) {
    stop("invalid argument(s)")
})


## extract the phylo4 part of phylo4d; relies on implicit coerce method


#' Get tree from tree+data object
#' 
#' Extracts a \code{phylo4} tree object from a \code{phylo4d} tree+data object.
#' 
#' \code{extractTree} extracts just the phylogeny from a tree+data object. The
#' phylogeny contains the topology (how the nodes are linked together), the
#' branch lengths (if any), and any tip and/or node labels. This may be useful
#' for extracting a tree from a \code{phylo4d} object, and associating with
#' another phenotypic dataset, or to convert the tree to another format.
#' 
#' @param from a \code{phylo4d} object, containing a phylogenetic tree plus
#' associated phenotypic data. Created by the \code{phylo4d()} function.
#' @author Ben Bolker
#' @seealso \code{\link{phylo4}}, \code{\link{phylo4d}},
#' \code{\link{coerce-methods}} for translation functions.
#' @keywords methods
#' @examples
#' 
#' tree.phylo <- ape::read.tree(text = "((a,b),c);")
#' tree <- as(tree.phylo, "phylo4")
#' plot(tree)
#' tip.data <- data.frame(size = c(1, 2, 3), row.names = c("a", "b", "c"))
#' (treedata <- phylo4d(tree, tip.data))
#' plot(treedata)
#' (tree1 <- extractTree(treedata))
#' plot(tree1)
#' 
extractTree <- function(from) {
    as(from, "phylo4")
}
