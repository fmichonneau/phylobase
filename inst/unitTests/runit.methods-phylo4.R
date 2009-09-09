#
# --- Test methods-phylo4.R ---
#
 
# create a phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,NA,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
elen <- descendant/10
elab <- paste("e", ancestor, descendant, sep="-")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# now alter internal ordering of each slot so nothing matches up;
# methods below should be able to handle this
phy@tip.label <- rev(phy@tip.label)
phy@node.label <- rev(phy@node.label)
phy@edge <- phy@edge[c(6:9, 1:5), ]
phy@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy@edge.label <- phy@edge.label[c(8:9, 1:7)]

# update test targets for edge-related slots
ancestor <- ancestor[c(6:9, 1:5)]
descendant <- descendant[c(6:9, 1:5)]
edge <- cbind(ancestor, descendant)
elen <- elen[c(6:9, 1:5)]
elab <- elab[c(6:9, 1:5)]

#-----------------------------------------------------------------------

test.nTips.phylo4 <- function() {
  checkEquals(nTips(phy), length(nid.tip))
}

test.nTips.ANY <- function() {
  # nTips phylo
  tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
  checkEquals(nTips(tr), 5)
}

test.nNodes.phylo4 <- function() {
  checkEquals(nNodes(phy), length(nid.int))
}

test.nodeType.phylo4 <- function() {
  checkIdentical(nodeType(phy), setNames(c(rep("tip", length(nid.tip)),
    "root", rep("internal", length(nid.int)-1)), c(nid.tip, nid.int)))
}

test.nodeId.phylo4 <- function() {
  checkIdentical(nodeId(phy), c(nid.tip, nid.int))
  checkIdentical(nodeId(phy, "all"), c(nid.tip, nid.int))
  checkIdentical(nodeId(phy, "tip"), nid.tip)
  checkIdentical(nodeId(phy, "internal"), nid.int)
  checkIdentical(nodeId(phy, "root"), nid.int[1])
}

test.nEdges.phylo4 <- function() {
  checkIdentical(nEdges(phy), nrow(edge))
}

test.edges.phylo4 <- function() {
  checkIdentical(edges(phy), edge)
  checkIdentical(edges(phy, drop.root=TRUE), edge[!is.na(edge[,1]),])
}

test.edgeOrder.phylo4 <- function() {
  checkIdentical(edgeOrder(phy), "unknown")
  checkIdentical(edgeOrder(reorder(phy, "preorder")), "preorder")
  checkIdentical(edgeOrder(reorder(phy, "postorder")), "postorder")
}

test.edgeId.phylo4 <- function() {
  eid <- paste(ancestor, descendant, sep="-")
  checkIdentical(edgeId(phy), eid)
  checkIdentical(edgeId(phy, "all"), eid)
  checkIdentical(edgeId(phy, "tip"), eid[descendant %in% nid.tip])
  checkIdentical(edgeId(phy, "internal"), eid[!descendant %in% nid.tip])
  checkIdentical(edgeId(phy, "root"), eid[is.na(ancestor)])
}

test.hasEdgeLength.phylo4 <- function() {
  checkTrue(hasEdgeLength(phy))
  phy@edge.length <- NA_real_
  checkTrue(!hasEdgeLength(phy))
}

test.edgeLength.phylo4 <- function() {
  # all edge lengths
  checkIdentical(edgeLength(phy), setNames(elen, paste(ancestor,
    descendant, sep="-")))
  # one edge length, by label
  checkEquals(edgeLength(phy, "t1"), c(`7-1`=0.1))
  # one edge length, by node ID
  checkEquals(edgeLength(phy, 1), c(`7-1`=0.1))
  # non-existent edge, by label
  ans <- structure(NA_real_, .Names = NA_character_)
  checkEquals(suppressWarnings(edgeLength(phy, "xxx")), ans)
  # non-existent edge, by number
  checkEquals(suppressWarnings(edgeLength(phy, 999)), ans)
}

test.Replace.edgeLength.phylo4 <- function() {
  #TODO function(x, use.names=TRUE, ..., value)
}

test.sumEdgeLength.phylo4 <- function() {
  #TODO function(phy, node)
}

test.isRooted.phylo4 <- function() {
  checkTrue(isRooted(phy))
}

test.rootNode.phylo4 <- function() {
  checkIdentical(rootNode(phy), nid.int[1])
}

test.Replace.rootNode.phylo4 <- function() {
  #TODO function(x, value)
}

test.labels.phylo4 <- function() {
  # function(object, type = c("all", "tip", "internal"))
  checkIdentical(labels(phy), setNames(c(lab.tip, lab.int), c(nid.tip,
    nid.int)))
  checkIdentical(labels(phy, "all"), setNames(c(lab.tip, lab.int),
    c(nid.tip, nid.int)))
  checkIdentical(labels(phy, "tip"), setNames(lab.tip, nid.tip))
  checkIdentical(labels(phy, "internal"), setNames(lab.int, nid.int))
}

test.Replace.labels.phylo4 <- function() {
  #TODO function(object, type = c("tip", "internal", "allnode"), use.names, ..., value)
}

test.hasNodeLabels.phylo4 <- function() {
  checkTrue(hasNodeLabels(phy))
  phy@node.label <- NA_character_
  checkTrue(!hasNodeLabels(phy))
}

test.nodeLabels.phylo4 <- function() {
  checkIdentical(nodeLabels(phy), setNames(lab.int, nid.int))
}

test.Replace.nodeLabels.phylo4 <- function() {
  #TODO function(object, ...,  value) {
}

test.tipLabels.phylo4 <- function() {
  checkIdentical(tipLabels(phy), setNames(lab.tip, nid.tip))
}

test.Replace.tipLabels.phylo4 <- function() {
  #TODO function(object, ...,  value) {
}

test.hasEdgeLabels.phylo4 <- function() {
  checkTrue(hasEdgeLabels(phy))
  phy@edge.label <- NA_character_
  checkTrue(!hasEdgeLabels(phy))
}

test.edgeLabels.phylo4 <- function() {
  checkIdentical(edgeLabels(phy), setNames(elab, paste(ancestor,
    descendant, sep="-")))
}

test.Replace.edgeLabels.phylo4 <- function() {
  #TODO function(object, ...,  value) {
}

test.print.phylo4 <- function() {
  #TODO? this just calls printphylo4 function
}

test.show.phylo4 <- function() {
  #TODO? this just calls printphylo4 function
}

test.names.phylo4 <- function() {
  #TODO?
}

test.head.phylo4 <- function() {
  #TODO?
}

test.tail.phylo4 <- function() {
  #TODO?
}

test.summary.phylo4 <- function() {
 #TODO? function (object, quiet=FALSE)
}

# not an exported function -- called internally by reorder("phylo4")
#test.orderIndex <- function() {
#}

test.reorder.phylo4 <- function() {
  #TODO
}


