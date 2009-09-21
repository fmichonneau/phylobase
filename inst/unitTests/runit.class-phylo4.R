#
# --- Test class-phylo4.R ---
#

test.phylo4.matrix <- function() {
    edge <- structure(c(6L, 7L, 8L, 8L, 9L, 9L, 7L, 6L, 7L, 8L, 1L, 9L,
      2L, 3L, 4L, 5L), .Dim = c(8, 2))
    edge.length <- c(0.2, 0.5, 0.2, 0.15, 0.1, 0.1, 0.7, 1)
    tip.label <- paste("t", 1:5, sep="")
    node.label <- paste("n", 1:4, sep="")
    edge.label <- paste("e", 1:8, sep="") 
    order="preorder"
    annote <- list(x="annotation")
    phy <- phylo4(edge, edge.length=edge.length, tip.label=tip.label,
      node.label=node.label, edge.label=edge.label, order=order,
      annote=annote)

    # test each slot
    checkIdentical(edge, unname(edges(phy)))
    checkIdentical(edge.length, unname(edgeLength(phy)))
    checkIdentical(4L, nNodes(phy))
    checkIdentical(tip.label, unname(tipLabels(phy)))
    checkIdentical(node.label, unname(nodeLabels(phy)))
    checkIdentical(edge.label, unname(edgeLabels(phy)))
    checkIdentical(order, edgeOrder(phy))
    checkIdentical(annote, phy@annote)

    # test improper cases
    #checkException(phylo4(edge, edge.length=999))  # recycling is allowed?
    checkException(phylo4(edge, tip.label=999))
    checkException(phylo4(edge, node.label=999))
    #checkException(phylo4(edge, edge.label=999))  # recycling is allowed?
    checkException(phylo4(edge, order="invalid order"))
    checkException(phylo4(edge, annote="invalid annotation"))
}

# note: this method mostly just wraps phylo->phylo4 coercion, which is
# tested more thoroughly in runit.setAs-methods.R; focus here is on
# annote and check.node.labels arguments
test.phylo4.phylo <- function() {
    tr <- read.tree(text="(((t1:0.2,(t2:0.1,t3:0.1):0.15):0.5,t4:0.7):0.2,t5:1):0.4;")

    ##
    ## annote
    ##

    annote <- list(x="annotation")
    phy <- phylo4(tr, annote=annote)
    checkIdentical(annote, phy@annote)

    ##
    ## check.node.labels
    ##

    # case 0: no node labels
    phy <- phylo4(tr)
    checkTrue(!hasNodeLabels(phy))

    # case 1: keep unique character labels
    tr$node.label <- paste("n", 1:4, sep="")
    phy <- phylo4(tr, check.node.labels="keep")
    checkIdentical(tr$node.label, unname(nodeLabels(phy)))
    # keeping node labels should be the default
    checkIdentical(phy, phylo4(tr))

    # case 2: keep unique number-like character labels
    tr$node.label <- as.character(1:4)
    phy <- phylo4(tr, check.node.labels="keep")
    checkIdentical(tr$node.label, unname(nodeLabels(phy)))

    # case 3: keep unique numeric labels, but convert to character
    tr$node.label <- as.numeric(1:4)
    phy <- phylo4(tr, check.node.labels="keep")
    checkIdentical(as.character(tr$node.label), unname(nodeLabels(phy)))

    # case 4: must drop non-unique labels
    tr$node.label <- rep("x", 4)
    checkException(phylo4(tr))
    checkException(phylo4(tr, check.node.labels="keep"))
    # test dropping node labels
    phy <- phylo4(tr, check.node.labels="drop")
    checkTrue(!hasNodeLabels(phy))

}
