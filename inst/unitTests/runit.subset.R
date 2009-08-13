#
# --- Test subset methods ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 

test.subset.phylo <- function() {
    print(subset(tr, 1:4))
    print(subset(tr, 1:2))
}

test.subset.phylo4 <- function() {
    phy <- as(tr, "phylo4")
    print(subset(phy, 1:4))
    print(subset(phy, 1:2))
}

test.subset.phylo4d <- function() {
    phyd <- as(tr, "phylo4d")
    print(subset(phyd, 1:4))
## the following print statement currently fails, for reasons related to
## the subtle differences between str(as(subset(tr, 1:2), "phylo4d"))
## and str(subset(as(tr, "phylo4d"), 1:2))
#    print(subset(phyd, 1:2))
}

