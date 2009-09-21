#
# --- Test treewalk.R ---
#

# Create sample phylo4 tree for testing
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
phy <- as(tr, "phylo4")

test.getNode <- function() {
# Note: we're not explicitly testing missing="warn" condition below;
# however, if "OK" and "fail" both work as expected, then so must "warn"
    # node only has valid characters
    checkEquals(getNode(phy, "spA"), c(spA=1))
    checkEquals(getNode(phy, c("spA", "spC")), c(spA=1, spC=3))

    # node only has valid integers 
    ans <- 4
    names(ans) <- "spD"
    checkEquals(getNode(phy, 4), ans)
    ans <- c(4,6)
    names(ans) <- c("spD", NA)
    checkEquals(getNode(phy, c(4,6)), ans)

    # node includes only missing characters (names), but missing=OK
    ans <- rep(NA_integer_, 2)  # return values should be NA
    names(ans) <- rep(NA, 2)  # return values should have NA names
    checkEquals(getNode(phy, c("xxx", "yyy"), missing="OK"), ans)
    # now missing = "fail"
    checkException(getNode(phy, c("xxx", "yyy"), missing="fail"))

    # node includes only missing numbers (IDs), but missing=OK
    ans <- rep(NA_integer_, 3)  # return values should be NA
    names(ans) <- rep(NA, 3)  # return values should have NA names
    checkEquals(getNode(phy, c(-9, 0, 50), missing="OK"), ans)
    # now missing = "fail"
    checkException(getNode(phy, c(-9, 0, 50), missing="fail"), ans)

    # node includes NAs, but missing = "OK"
    checkTrue(is.na(getNode(phy, NA_integer_, missing="OK")))
    checkTrue(is.na(getNode(phy, NA_character_, missing="OK")))

    # node includes mixture of valid values and NAs
    ans <- c(2, NA)
    names(ans) <- c("spB", NA) 
    checkEquals(getNode(phy, c("spB", NA), missing="OK"), ans)
    checkEquals(getNode(phy, c(2, NA), missing="OK"), ans)

    # node is neither integer-like nor character
    checkException(getNode(phy, 1.5))
}

test.ancestor <- function() {
    # function(phy,node)
}

test.children <- function() {
    # function(phy,node)
}

test.descendants <- function() {
    # function (phy, node, type=c("tips","children","all"))
    phy <- phylo4(read.tree(text="((t3,t4),(t1,(t2,t5)));"))

    # node = tip
    checkIdentical(descendants(phy, 5),
        setNames(5L, "t5"))
    checkIdentical(descendants(phy, 5, "tips"),
        setNames(5L, "t5"))
    checkIdentical(descendants(phy, 5, "children"),
        setNames(integer(0), character(0)))
    checkIdentical(descendants(phy, 5, "all"),
        setNames(5L, "t5"))

    # node = internal
    checkIdentical(descendants(phy, 8),
        setNames(c(3L, 4L, 5L), c("t1", "t2", "t5")))
    checkIdentical(descendants(phy, 8, "tips"),
        setNames(c(3L, 4L, 5L), c("t1", "t2", "t5")))
    checkIdentical(descendants(phy, 8, "children"),
        setNames(c(3L, 9L), c("t1", NA)))
    checkIdentical(descendants(phy, 8, "all"),
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
}


