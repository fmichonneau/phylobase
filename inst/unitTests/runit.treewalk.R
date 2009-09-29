#
# --- Test treewalk.R ---
#

# Create sample phylo4 tree for testing
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
phytr <- as(tr, "phylo4")

# create phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
nid.all <- c(nid.tip, nid.int)
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
lab.all <- c(lab.tip, lab.int)
eid <- paste(ancestor, descendant, sep="-")
elen <- descendant/10
elab <- paste("e", eid, sep="")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@label <- rev(phy@label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

# update test targets for edge-related slots
ancestor <- ancestor[c(6:9, 1:5)]
descendant <- descendant[c(6:9, 1:5)]
edge <- cbind(ancestor, descendant)
eid <- eid[c(6:9, 1:5)]
elen <- elen[c(6:9, 1:5)]
elab <- elab[c(6:9, 1:5)]

#-----------------------------------------------------------------------

test.getNode <- function() {
# Note: we're not explicitly testing missing="warn" condition below;
# however, if "OK" and "fail" both work as expected, then so must "warn"
    # node only has valid characters
    checkEquals(getNode(phytr, "spA"), c(spA=1))
    checkEquals(getNode(phytr, c("spA", "spC")), c(spA=1, spC=3))

    # node only has valid integers 
    ans <- 4
    names(ans) <- "spD"
    checkEquals(getNode(phytr, 4), ans)
    ans <- c(4,6)
    names(ans) <- c("spD", NA)
    checkEquals(getNode(phytr, c(4,6)), ans)

    # node includes only missing characters (names), but missing=OK
    ans <- rep(NA_integer_, 2)  # return values should be NA
    names(ans) <- rep(NA, 2)  # return values should have NA names
    checkEquals(getNode(phytr, c("xxx", "yyy"), missing="OK"), ans)
    # now missing = "fail"
    checkException(getNode(phytr, c("xxx", "yyy"), missing="fail"))

    # node includes only missing numbers (IDs), but missing=OK
    ans <- rep(NA_integer_, 3)  # return values should be NA
    names(ans) <- rep(NA, 3)  # return values should have NA names
    checkEquals(getNode(phytr, c(-9, 0, 50), missing="OK"), ans)
    # now missing = "fail"
    checkException(getNode(phytr, c(-9, 0, 50), missing="fail"), ans)

    # node includes NAs, but missing = "OK"
    checkTrue(is.na(getNode(phytr, NA_integer_, missing="OK")))
    checkTrue(is.na(getNode(phytr, NA_character_, missing="OK")))

    # node includes mixture of valid values and NAs
    ans <- c(2, NA)
    names(ans) <- c("spB", NA) 
    checkEquals(getNode(phytr, c("spB", NA), missing="OK"), ans)
    checkEquals(getNode(phytr, c(2, NA), missing="OK"), ans)

    # node is neither integer-like nor character
    checkException(getNode(phytr, 1.5))
}

test.ancestor <- function() {
    # function(phy,node)
}

test.children <- function() {
    # function(phy,node)
}

test.descendants <- function() {
    # function (phy, node, type=c("tips","children","all"))
    phytr <- phylo4(read.tree(text="((t3,t4),(t1,(t2,t5)));"))

    # node = tip
    checkIdentical(descendants(phytr, 5),
        setNames(5L, "t5"))
    checkIdentical(descendants(phytr, 5, "tips"),
        setNames(5L, "t5"))
    checkIdentical(descendants(phytr, 5, "children"),
        setNames(integer(0), character(0)))
    checkIdentical(descendants(phytr, 5, "all"),
        setNames(5L, "t5"))

    # node = internal
    checkIdentical(descendants(phytr, 8),
        setNames(c(3L, 4L, 5L), c("t1", "t2", "t5")))
    checkIdentical(descendants(phytr, 8, "tips"),
        setNames(c(3L, 4L, 5L), c("t1", "t2", "t5")))
    checkIdentical(descendants(phytr, 8, "children"),
        setNames(c(3L, 9L), c("t1", NA)))
    checkIdentical(descendants(phytr, 8, "all"),
        setNames(c(3L, 9L, 4L, 5L), c("t1", NA, "t2", "t5")))
}

test.siblings <- function() {
    # function(phy, node, include.self=FALSE)
}

test.ancestors <- function() {
    # function (phy, node, type=c("all","parent","ALL"))
}

test.MRCA <- function() {
    # function(phy, ...)
}

test.shortestPath <- function() {
    # function(phy, node1, node2)
}

test.getEdge <- function() {
    # function(phy, node, type=c("descendant", "ancestor"),
    #     missing=c("warn", "OK", "fail"))

    #
    # nodes as descendants
    #

    # node only has valid descendants, as characters
    checkIdentical(getEdge(phy.alt, "t1"), setNames("7-1", 1))
    checkIdentical(getEdge(phy.alt, c("t1", "t3")), setNames(c("7-1",
        "8-3"), c(1,3)))

    # node only has valid descendants, as integers 
    checkIdentical(getEdge(phy.alt, 1), setNames("7-1", 1))
    checkIdentical(getEdge(phy.alt, c(1,3)), setNames(c("7-1",
        "8-3"), c(1,3)))

    # node includes only missing characters (labels), but missing=OK
    checkIdentical(getEdge(phy.alt, c("x", "y", "z"), missing="OK"),
        setNames(rep(NA, 3), rep(NA, 3)))
    # now missing = "fail"
    checkException(getEdge(phy.alt, c("x", "y", "z"), missing="fail"))

    # node includes only missing numbers (IDs), but missing=OK
    checkIdentical(getEdge(phy.alt, c(-9, 0, 50), missing="OK"),
        setNames(rep(NA, 3), rep(NA, 3)))
    # now missing = "fail"
    checkException(getEdge(phy, c(-9, 0, 50), missing="fail"))

    # node includes NAs, but missing = "OK"
    checkTrue(is.na(getEdge(phy, NA_integer_, missing="OK")))
    checkTrue(is.na(getEdge(phy, NA_character_, missing="OK")))

    # node includes mixture of valid values and NAs
    checkIdentical(getEdge(phy, c("t3", NA), missing="OK"),
        setNames(c("8-3", NA), c(3, NA)))
    checkIdentical(getEdge(phy, c(3, NA), missing="OK"),
        setNames(c("8-3", NA), c(3, NA)))

    # node is neither integer-like nor character
    checkException(getEdge(phy, 1.5))

    #
    # nodes as ancestors
    #

    # node only has valid ancestors, as characters
    checkIdentical(getEdge(phy.alt, "n6", type="ancestor"),
        setNames(c("6-7", "6-8"), c(6, 6)))
    checkIdentical(getEdge(phy.alt, c("n6", "n8"), type="ancestor"),
        setNames(c("6-7", "6-8", "8-9", "8-3"), c(6, 6, 8, 8)))

    # node only has valid ancestors, as integers 
    checkIdentical(getEdge(phy.alt, 6, type="ancestor"),
        setNames(c("6-7", "6-8"), c(6, 6)))
    checkIdentical(getEdge(phy.alt, c(6, 8), type="ancestor"),
        setNames(c("6-7", "6-8", "8-9", "8-3"), c(6, 6, 8, 8)))

    # node includes only missing characters (labels), but missing=OK
    checkIdentical(getEdge(phy.alt, c("x", "y", "z"), type="ancestor",
        missing="OK"), setNames(rep(NA, 3), rep(NA, 3)))
    # node includes only tips (labels), but missing=OK
    checkIdentical(getEdge(phy.alt, c("t1", "t3"), type="ancestor",
        missing="OK"), setNames(rep(NA, 2), c(1, 3)))
    # now missing = "fail"
    checkException(getEdge(phy.alt, c("x", "y", "z"), missing="fail"))
    checkException(getEdge(phy.alt, c("t1", "t3"), type="ancestor",
        missing="fail"))

    # node includes only missing numbers (IDs), but missing=OK
    checkIdentical(getEdge(phy.alt, c(-9, 0, 50), type="ancestor",
        missing="OK"), setNames(rep(NA, 3), rep(NA, 3)))
    # node includes only tips (labels), but missing=OK
    checkIdentical(getEdge(phy.alt, c(1, 3), type="ancestor",
        missing="OK"), setNames(rep(NA, 2), c(1, 3)))
    # now missing = "fail"
    checkException(getEdge(phy.alt, c(-9, 0, 50), missing="fail"))
    checkException(getEdge(phy.alt, c(1, 3), type="ancestor",
        missing="fail"))

    # node includes NAs, but missing = "OK"
    checkTrue(is.na(getEdge(phy.alt, NA_integer_, type="ancestor",
        missing="OK")))
    checkTrue(is.na(getEdge(phy.alt, NA_character_, type="ancestor",
        missing="OK")))

    # node includes mixture of valid values and NAs
    checkIdentical(getEdge(phy.alt, c("t3", "n8", NA), type="ancestor",
        missing="OK"), setNames(c(NA, "8-9", "8-3", NA), c(3, 8, 8, NA)))
    checkIdentical(getEdge(phy.alt, c(3, 8, NA), type="ancestor",
        missing="OK"), setNames(c(NA, "8-9", "8-3", NA), c(3, 8, 8, NA)))

}


