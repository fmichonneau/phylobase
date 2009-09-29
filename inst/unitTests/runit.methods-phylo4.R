#
# --- Test methods-phylo4.R ---
#
 
# create ape::phylo version of a simple tree for testing
nwk <- "((t1:0.1,t2:0.2)n7:0.7,(t3:0.3,(t4:0.4,t5:0.5)n9:0.9)n8:0.8)n6:0.6;"
tr <- read.tree(text=nwk)

# create analogous phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
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

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@tip.label <- rev(phy@tip.label)
phy.alt@node.label <- rev(phy@node.label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

# update test targets for edge-related slots
ancestor <- ancestor[c(6:9, 1:5)]
descendant <- descendant[c(6:9, 1:5)]
edge <- cbind(ancestor, descendant)
elen <- elen[c(6:9, 1:5)]
elab <- elab[c(6:9, 1:5)]

#-----------------------------------------------------------------------

test.nTips.phylo4 <- function() {
  checkEquals(nTips(phy.alt), length(nid.tip))
}

test.nTips.ANY <- function() {
  # nTips phylo
  checkEquals(nTips(tr), 5)
}

test.nNodes.phylo4 <- function() {
  checkEquals(nNodes(phy.alt), length(nid.int))
}

test.nodeType.phylo4 <- function() {
  checkIdentical(nodeType(phy.alt), setNames(c(rep("tip", length(nid.tip)),
    "root", rep("internal", length(nid.int)-1)), c(nid.tip, nid.int)))
}

test.nodeId.phylo4 <- function() {
  checkIdentical(nodeId(phy.alt), c(nid.tip, nid.int))
  checkIdentical(nodeId(phy.alt, "all"), c(nid.tip, nid.int))
  checkIdentical(nodeId(phy.alt, "tip"), nid.tip)
  checkIdentical(nodeId(phy.alt, "internal"), nid.int)
  checkIdentical(nodeId(phy.alt, "root"), nid.int[1])
}

test.nEdges.phylo4 <- function() {
  checkIdentical(nEdges(phy.alt), nrow(edge))
}

test.edges.phylo4 <- function() {
  checkIdentical(edges(phy.alt), edge)
  checkIdentical(edges(phy.alt, drop.root=TRUE), edge[edge[,1] != 0,])
}

test.edgeOrder.phylo4 <- function() {
  checkIdentical(edgeOrder(phy.alt), "unknown")
  checkIdentical(edgeOrder(reorder(phy.alt, "preorder")), "preorder")
  checkIdentical(edgeOrder(reorder(phy.alt, "postorder")), "postorder")
}

test.edgeId.phylo4 <- function() {
  eid <- paste(ancestor, descendant, sep="-")
  checkIdentical(edgeId(phy.alt), eid)
  checkIdentical(edgeId(phy.alt, "all"), eid)
  checkIdentical(edgeId(phy.alt, "tip"), eid[descendant %in% nid.tip])
  checkIdentical(edgeId(phy.alt, "internal"), eid[!descendant %in% nid.tip])
  checkIdentical(edgeId(phy.alt, "root"), eid[ancestor == 0])
}

test.hasEdgeLength.phylo4 <- function() {
  checkTrue(hasEdgeLength(phy.alt))
  phy.alt@edge.length <- NA_real_
  checkTrue(!hasEdgeLength(phy.alt))
}

test.edgeLength.phylo4 <- function() {
  # all edge lengths
  checkIdentical(edgeLength(phy.alt), setNames(elen, paste(ancestor,
    descendant, sep="-")))
  # one edge length, by label
  checkEquals(edgeLength(phy.alt, "t1"), c(`7-1`=0.1))
  # one edge length, by node ID
  checkEquals(edgeLength(phy.alt, 1), c(`7-1`=0.1))
  # non-existent edge, by label
  ans <- structure(NA_real_, .Names = NA_character_)
  checkEquals(suppressWarnings(edgeLength(phy.alt, "xxx")), ans)
  # non-existent edge, by number
  checkEquals(suppressWarnings(edgeLength(phy.alt, 999)), ans)
}

test.Replace.edgeLength.phylo4 <- function() {
  #TODO function(x, use.names=TRUE, ..., value)
}

test.sumEdgeLength.phylo4 <- function() {
  #TODO function(phy, node)
}

test.isRooted.phylo4 <- function() {
  checkTrue(isRooted(phy.alt))
}

test.rootNode.phylo4 <- function() {
  checkIdentical(rootNode(phy.alt), nid.int[1])
}

test.Replace.rootNode.phylo4 <- function() {
  #TODO function(x, value)
}

test.labels.phylo4 <- function() {
  # function(object, type = c("all", "tip", "internal"))
  checkIdentical(labels(phy.alt), setNames(c(lab.tip, lab.int), c(nid.tip,
    nid.int)))
  checkIdentical(labels(phy.alt, "all"), setNames(c(lab.tip, lab.int),
    c(nid.tip, nid.int)))
  checkIdentical(labels(phy.alt, "tip"), setNames(lab.tip, nid.tip))
  checkIdentical(labels(phy.alt, "internal"), setNames(lab.int, nid.int))
}

test.Replace.labels.phylo4 <- function() {
  #TODO function(object, type = c("tip", "internal", "allnode"), use.names, ..., value)
}

test.hasNodeLabels.phylo4 <- function() {
  checkTrue(hasNodeLabels(phy.alt))
  phy.alt@node.label <- NA_character_
  checkTrue(!hasNodeLabels(phy.alt))
}

test.nodeLabels.phylo4 <- function() {
  checkIdentical(nodeLabels(phy.alt), setNames(lab.int, nid.int))
}

test.Replace.nodeLabels.phylo4 <- function() {
  #TODO function(object, ...,  value) {
}

test.tipLabels.phylo4 <- function() {
  checkIdentical(tipLabels(phy.alt), setNames(lab.tip, nid.tip))
}

test.Replace.tipLabels.phylo4 <- function() {
  #TODO function(object, ...,  value) {
}

test.hasEdgeLabels.phylo4 <- function() {
  checkTrue(hasEdgeLabels(phy.alt))
  phy.alt@edge.label <- NA_character_
  checkTrue(!hasEdgeLabels(phy.alt))
}

test.edgeLabels.phylo4 <- function() {
  checkIdentical(edgeLabels(phy.alt), setNames(elab, paste(ancestor,
    descendant, sep="-")))
}

test.Replace.edgeLabels.phylo4 <- function() {
  #TODO function(object, ...,  value) {
}

## this is also the print method
## this mostly just wraps .phylo4ToDataFrame, which is tested elsewhere
##test.show.phylo4 <- function() {
##}

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
  phy.sum <- summary(phy.alt, quiet=TRUE)
  checkIdentical(phy.sum$name, "phy.alt")
  checkIdentical(phy.sum$nb.tips, length(nid.tip))
  checkIdentical(phy.sum$nb.nodes, length(nid.int))
  checkIdentical(phy.sum$mean.el, mean(elen))
  checkIdentical(phy.sum$var.el, var(elen))
  checkIdentical(phy.sum$sumry.el, summary(elen))
  # now make root edge length NA
  edgeLength(phy.alt)[edgeId(phy.alt, "root")] <- NA
  phy.sum2 <- summary(phy.alt, quiet=TRUE)
  checkIdentical(phy.sum2$mean.el, mean(edgeLength(phy.alt), na.rm=TRUE))
  checkIdentical(phy.sum2$var.el, var(edgeLength(phy.alt), na.rm=TRUE))
  checkIdentical(phy.sum2$sumry.el, summary(na.omit(edgeLength(phy.alt))))
  # now remove edge lengths altogether
  phy.alt@edge.length[] <- NA
  phy.sum3 <- summary(phy.alt, quiet=TRUE)
  checkTrue(is.null(phy.sum3$mean.el))
  checkTrue(is.null(phy.sum3$var.el))
  checkTrue(is.null(phy.sum3$sumry.el))
}

# not an exported function -- called internally by reorder("phylo4")
#test.orderIndex <- function() {
#}

test.reorder.phylo4 <- function() {
  #TODO
}


