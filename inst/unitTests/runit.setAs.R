#
# --- Test setAs methods ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
phy <- as(tr, "phylo4")

test.phylo4.As.data.frame <- function() {

    # rooted tree
    checkTrue(is.data.frame(as(phy, "data.frame")))
    phy.df <- structure(list(label = c(NA, NA, NA, NA, "spA", "spB",
       "spC", "spD", "spE"), node = c(6L, 7L, 8L, 9L, 1L, 2L, 3L, 4L,
        5L), ancestor = c(NA, 6L, 7L, 8L, 8L, 9L, 9L, 7L, 6L),
        edge.length = c(0.4, 0.2, 0.5, 0.15, 0.2, 0.1, 0.1, 0.7, 1),
        node.type = structure(c(2L, 1L, 1L, 1L, 3L, 3L, 3L, 3L, 3L),
        .Label = c("internal", "root", "tip"), class = "factor")),
        .Names = c("label", "node", "ancestor", "edge.length",
        "node.type"), row.names = c(6L, 7L, 8L, 9L, 1L, 2L, 3L, 4L, 5L),
        class = "data.frame")
    checkEquals(as(phy, "data.frame"), phy.df)

    # unrooted tree
    tru <- unroot(tr)
    phyu <- as(tru, "phylo4")
    # should probably check that this coercion results in something
    # *correct*, not just that it produces a data.frame
    checkTrue(is.data.frame(as(phyu, "data.frame")))

}
