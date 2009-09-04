#
# --- Test prune.R ---
#

data(geospiza)
gtree <- extractTree(geospiza)

test.prune.phylo4 <- function() {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
    checkIdentical(gtree, prune(gtree, character(0)))
}

test.prune.phylo4d <- function() {
    # function(phy, tip, trim.internal = TRUE, subtree = FALSE, ...)
    checkIdentical(geospiza, prune(geospiza, character(0)))
}
