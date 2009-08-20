#
# --- Test phylo4 methods ---
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

test.nTips.phylo4 <- function() {
  checkEquals(5, nTips(phy))
  checkEquals(5, nTips(phyd))
}

test.nTips.ANY <- function() {
  # nTips phylo
  checkEquals(5, nTips(tr))
}

test.nNodes.phylo4 <- function() {
  checkEquals(4, nNodes(phy))
  checkEquals(4, nNodes(phyd))
}

test.nodeType.phylo4 <- function() {
}

test.nodeId.phylo4 <- function() {
  # do for type=c("internal","tip","allnode")
}

test.nEdges.phylo4 <- function() {
}

test.edges.phylo4 <- function() {
  # function(x, order, drop.root=FALSE, ...)
}

test.edgeOrder.phylo4 <- function() {
  # function(x, ...)
}

test.hasEdgeLength.phylo4 <- function() {
  checkTrue(hasEdgeLength(phy))
}

test.edgeLength.phylo4 <- function() {
  # all edge lengths
  ans <- structure(c(0.2, 0.5, 0.2, 0.15, 0.1, 0.4, 0.1, 0.7, 1),
    .Names = c("6-7", "7-8", "8-1", "8-9", "9-2", "NA-6", "9-3", "7-4",
    "6-5"))
  checkEquals(ans, edgeLength(phy))
  # one edge length, by label
  ans <- c(`8-1`=0.2)
  checkEquals(ans, edgeLength(phy, "spA"))
  # one edge length, by number
  checkEquals(ans, edgeLength(phy, 1))
  # non-existent edge, by label
  ans <- structure(NA_real_, .Names = NA_character_)
  checkEquals(ans, edgeLength(phy, "xxx"))
  # non-existent edge, by number
  checkEquals(ans, edgeLength(phy, 999))
}

test.Replace.edgeLength.phylo4 <- function() {
  # function(x, use.names=TRUE, ..., value)
}

test.sumEdgeLength.phylo4 <- function() {
  # function(phy, node)
}

test.isRooted.phylo4 <- function() {
}

test.rootNode.phylo4 <- function() {
}

test.Replace.rootNode.phylo4 <- function() {
  # function(x, value)
}

test.labels.phylo4 <- function() {
  # function(object, type = c("tip", "internal", "allnode"), ...)
}

test.Replace.labels.phylo4 <- function() {
  # signature(object="phylo4", type="ANY", use.names="ANY", value="character"),
  # function(object, type = c("tip", "internal", "allnode"), use.names, ..., value)
}

test.hasNodeLabels.phylo4 <- function() {
}

test.nodeLabels.phylo4 <- function() {
}

test.Replace.nodeLabels.phylo4 <- function() {
  # signature(object="phylo4", value="character")
  # function(object, ...,  value) {
}

test.tipLabels.phylo4 <- function() {
}

test.Replace.tipLabels.phylo4 <- function() {
  # signature(object="phylo4", value="character")
  # function(object, ...,  value) {
}

test.hasEdgeLabels.phylo4 <- function() {
}

test.edgeLabels.phylo4 <- function() {
}

test.Replace.edgeLabels.phylo4 <- function() {
  # signature(object="phylo4", value="character")
  # function(object, ...,  value) {
}

test.print.phylo4 <- function() {
  # this just calls printphylo4 function
}

test.show.phylo4 <- function() {
  # this just calls printphylo4 function
}

test.names.phylo4 <- function() {
}

test.head.phylo4 <- function() {
}

test.tail.phylo4 <- function() {
}

test.summary.phylo4 <- function() {
 # function (object, quiet=FALSE)
}

test.orderIndex <- function() {
  # function(phy, order = c('preorder', 'postorder'))
}

test.reorder.phylo4 <- function() {
}


