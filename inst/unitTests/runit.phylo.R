#
# --- Test ape import and handling ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 

test.phylo.to.phylo4.simple <- function() {
  phy <- as(tr, "phylo4")
  checkTrue(class(phy)=="phylo4")
  checkEquals(phy, phylo4(tr))
}

test.phylo.to.phylo4d.simple <- function() {
  phyd <- as(tr, "phylo4d")
  checkTrue(class(phyd)=="phylo4d")
  checkEquals(phyd, phylo4d(tr))
}

test.roundtrip.phylo.to.phylo4 <- function() {
  phy <- as(tr, "phylo4")
  checkEquals(tr, as(phy, "phylo"))
}

test.roundtrip.phylo.to.phylo4d <- function() {
  phyd <- as(tr, "phylo4d")
  checkEquals(tr, as(phyd, "phylo"))
}

test.phylo.import.with.character.node.labels <- function() {

  # case 1: unique non-numeric characters
  tr$node.label <- paste("n", 1:4, sep="")

    # import to phylo4
    tmp <- phylo4(tr, check.node.labels="keep")
    checkEquals(tmp, phylo4(tr))

    # import to phylo4d
    tmp <- phylo4d(tr, check.node.labels="keep")
    checkEquals(tmp, phylo4d(tr))
    checkEquals(unname(nodeLabels(tmp)), tr$node.label)
    checkEquals(nrow(tdata(tmp)), 0)

  # case 2: unique number-like characters
  tr$node.label <- as.character(1:4)

    # import to phylo4
    tmp <- phylo4(tr, check.node.labels="keep")
    checkEquals(tmp, phylo4(tr))

    # import to phylo4d
    tmp <- phylo4d(tr, check.node.labels="keep")
    checkEquals(tmp, phylo4d(tr))
    checkEquals(unname(nodeLabels(tmp)), tr$node.label)
    checkEquals(nrow(tdata(tmp)), 0)

  # case 3: non-unique characters
  tr$node.label <- rep("x", 4)

    # import to phylo4
    checkException(phylo4(tr))
    checkException(phylo4(tr, check.node.labels="keep"))

    # import to phylo4d
    checkException(phylo4d(tr))
    checkException(phylo4d(tr, check.node.labels="keep"))

}

test.phylo.import.with.numeric.node.labels <- function() {

    tr$node.label <- 1:4

    # keeping node labels should be the default
    checkEquals(phylo4(tr), phylo4(tr, check.node.labels="keep"))
    checkEquals(phylo4d(tr), phylo4d(tr, check.node.labels="keep"))

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

