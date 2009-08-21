#
# --- Test methods-phylo4.R ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
phy <- as(tr, "phylo4")
##   label node ancestr edge.length node.type
## 6  <NA>    6      NA        0.40      root
## 7  <NA>    7       6        0.20  internal
## 8  <NA>    8       7        0.50  internal
## 9  <NA>    9       8        0.15  internal
## 1   spA    1       8        0.20       tip
## 2   spB    2       9        0.10       tip
## 3   spC    3       9        0.10       tip
## 4   spD    4       7        0.70       tip
## 5   spE    5       6        1.00       tip
phyd <- as(phy, "phylo4d")

test.tdata.phylo4d <- function() {
    # function(x, type=c("tip", "internal", "allnode"),
    #   label.type=c("row.names","column"), empty.columns=TRUE, ...)
}

test.Replace.tdata.phylo4d <- function() {
    # function(object, type = c("tip", "internal", "allnode"), ...,
    #   value)
}

test.addData.phylo4d <- function() {
    # function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
    #   pos=c("after", "before"), merge.data=TRUE, match.data=TRUE, ...)
}

test.addData.phylo4 <- function() {
    # function(x, tip.data=NULL, node.data=NULL, all.data=NULL,
    #   pos=c("after", "before"), merge.data=TRUE, match.data=TRUE, ...)
}

test.summary.phylo4d <- function() {
}

test.hasNodeData.phylo4d <- function() {
}

test.na.omit.phylo4d <- function() {
    # function(object, ...)
}

