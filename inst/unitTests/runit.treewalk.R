#
# --- Test treewalk functions ---
#

# Note: we're not explicitly testing missing="warn" condition below;
# however, if "OK" and "fail" both work as expected, then so must "warn"
 
# Create sample phylo4 tree for testing
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
phy <- as(tr, "phylo4")

test.getNode.valid.character <- function() {
    # node only has valid characters
    checkEquals(getNode(phy, "spA"), c(spA=1))
    checkEquals(getNode(phy, c("spA", "spC")), c(spA=1, spC=3))
}

test.getNode.valid.integer <- function() {
    # node only has valid integers 
    ans <- 4
    names(ans) <- "spD"
    checkEquals(getNode(phy, 4), ans)
    ans <- c(4,6)
    names(ans) <- c("spD", NA)
    checkEquals(getNode(phy, c(4,6)), ans)
}

test.getNode.missing.character <- function() {
    # node includes only missing characters (names), but missing=OK
    ans <- rep(NA_integer_, 2)  # return values should be NA
    names(ans) <- c("xxx", "yyy")
    checkEquals(getNode(phy, c("xxx", "yyy"), missing="OK"), ans)
    # now missing = "fail"
    checkException(getNode(phy, c("xxx", "yyy"), missing="fail"))
}

test.getNode.missing.integer <- function() {
    # node includes only missing numbers (IDs), but missing=OK
    ans <- rep(NA_integer_, 3)  # return values should be NA
    names(ans) <- rep(NA, 3)
    checkEquals(getNode(phy, c(-9, 0, 50), missing="OK"), ans)
    # now missing = "fail"
    checkException(getNode(phy, c(-9, 0, 50), missing="fail"), ans)
}

test.getNode.NAs <- function() {
    # node includes NAs, but missing = "OK"
    checkTrue(is.na(getNode(phy, NA_integer_, missing="OK")))
    checkTrue(is.na(getNode(phy, NA_character_, missing="OK")))
}

test.getNode.mixed.cases <- function() {
    # node includes mixture of valid values and NAs
    ans <- c(2, NA)
    names(ans) <- c("spB", NA) 
    checkEquals(getNode(phy, c("spB", NA), missing="OK"), ans)
    checkEquals(getNode(phy, c(2, NA), missing="OK"), ans)
}

test.getNode.invalid.nodes <- function() {
    # node is neither integer-like nor character
    checkException(getNode(phy, 1.5))
}

