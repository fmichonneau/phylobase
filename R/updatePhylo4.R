updatePhylo4 <- function(phy, ...) {
    ## Add internal names for tip labels
    if(is.null(names(phy@tip.label))) {
        if(length(phy@tip.label == nTips(phy))) {
            names(phy@tip.label) <- nodeId(phy, "tip")
        }
        else stop("You have a problem with your tip labels")
    }

    ## Add internal names for node labels
    if(is.null(names(phy@node.label))) {
        if(length(phy@node.label) == nNodes(phy)) {
            names(phy@node.label) <- nodeId(phy, "internal")
        }
        else stop("You have a problem with your node labels.")
    }

    ## Add internal names for edge lengths
    if(hasEdgeLength(phy) && is.null(names(phy@edge.length))) {
        names(phy@edge.length) <- paste(phy@edge[,1], phy@edge[,2], sep="-")
    }

    if(is.character(msg <- checkPhylo4(phy))) stop(msg)
    else return(phy)

}
