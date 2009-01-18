## Alternative print method for phylo4, showing the contents of the tree data.
##  Not sure if it works for unrooted trees
printphylo <- function (x,printlen=6,...) {
    .Deprecated("print", package="phylobase")
    printlen <- max(1,printlen)
    nb.tip <- length(x@tip.label)
    nb.node <- x@Nnode
    nb.edge <- length(x@edge.label)
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and",
              nb.node, "internal nodes\n"))

    ## print tip labels
    cat("\nTips labels:\n")
    if (nb.tip > printlen) {
        cat(paste("\t", paste(x@tip.label[1:printlen], collapse = ", "),
                  ", ...\n", sep = ""))
    } else print(x@tip.label)

    ## print node labels
    cat("\nNodes labels:\n")
    if (nb.node > printlen) {
        cat(paste("\t", paste(x@node.label[1:printlen], collapse = ", "),
                  ", ...\n", sep = ""))
    } else print(x@node.label)

    ## print edge labels
    cat("\nEdges labels:\n")
    if (nb.edge > printlen) {
        cat(paste("\t", paste(x@edge.label[1:printlen], collapse = ", "),
                  ", ...\n", sep = ""))
    } else print(x@edge.label)

    ## slots
    ##     cat("\nSlots:\n")
    ##     cat(paste("@", names(x)[1:4], sep=""),sep="\t")
    ##     cat("\n")
    ##     cat(paste("@", names(x)[5:7], sep=""),sep="\t")
    ##     cat("\n")

    rlab <- if (isRooted(x)) "Rooted"  else "Unrooted"
    cat("\n", rlab, "; ", sep = "")
    blen <- if (hasEdgeLength(x))
      "includes branch lengths"
    else       "no branch lengths"
    cat(blen, "\n\n", sep = "")
}
