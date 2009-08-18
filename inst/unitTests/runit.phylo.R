#
# --- Test ape import and handling ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 

test.phylo.import.simple <- function() {
    checkTrue(class(phylo4(tr))=="phylo4")
    checkTrue(class(phylo4d(tr))=="phylo4d")
}

test.phylo.import.with.valid.node.labels <- function() {

    tr$node.label <- as.character(1:4)

    # import to phylo4
    tmp <- phylo4(tr, check.node.labels="keep")
    checkEquals(tmp, phylo4(tr))

    # import to phylo4d
    tmp <- phylo4d(tr, check.node.labels="keep")
    checkEquals(tmp, phylo4d(tr))
    checkEqualsNumeric(tmp@node.label, as.character(1:4))
    checkEquals(nrow(tdata(tmp)), 0)

}

test.phylo.import.with.numeric.node.labels <- function() {

    tr$node.label <- 1:4

    # can't keep invalid node labels
# TODO: remove/modify these after fm-branch merge
#    checkException(phylo4(tr))
#    checkException(phylo4(tr, check.node.labels="keep"))
#    checkException(phylo4d(tr))
#    checkException(phylo4d(tr, check.node.labels="keep"))

    # import to phylo4, dropping node labels
    tmp <- phylo4(tr, check.node.labels="drop")
    checkTrue(all(is.na(tmp@node.label)))

    # import to phylo4d, dropping node labels
    tmp <- phylo4d(tr, check.node.labels="drop")
    checkTrue(all(is.na(tmp@node.label)))
    checkEquals(nrow(tdata(tmp)), 0)

    # import to phylo4d, converting node labels to data
    tmp <- phylo4d(tr, check.node.labels="asdata")
    checkEquals(tdata(tmp, "internal", label.type="column")$labelValues,
        1:4)
    checkTrue(all(is.na(tmp@node.label)))
}

test.phylo.import.2tips <- function() {
    tr2 <- drop.tip(tr, 3:Ntip(tr))
    phy2 <- as(tr2, "phylo4")
    checkEquals(nTips(as(tr2, "phylo4")), 2)
    checkEquals(nNodes(as(tr2, "phylo4")), 1)
}

