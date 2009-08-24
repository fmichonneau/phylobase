#
# --- Test prune.R ---
#

data(geospiza)
gtree <- extractTree(geospiza)


test.DropTip <- function() {
    # function(phy, tip, ...)
}

test.prune.phylo4 <- function() {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
    checkEquals(gtree, prune(gtree, character(0)))
}

test.prune.phylo4d <- function() {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
    checkEquals(geospiza, prune(geospiza, character(0)))
}

test.prune.phylo <- function() {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
}
