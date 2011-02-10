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
nid.all <- c(nid.tip, nid.int)
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
lab.all <- c(lab.tip, lab.int)
eid <- paste(ancestor, descendant, sep="-")
elen <- descendant/10
elab <- paste("e", eid, sep="")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@label <- rev(phy@label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

# update test targets for edge-related slots
ancestor <- ancestor[c(6:9, 1:5)]
descendant <- descendant[c(6:9, 1:5)]
edge <- cbind(ancestor, descendant)
eid <- eid[c(6:9, 1:5)]
elen <- elen[c(6:9, 1:5)]
elab <- elab[c(6:9, 1:5)]

op <- phylobase.options()
#-----------------------------------------------------------------------

test.nTips.phylo4 <- function() {
  checkEquals(nTips(phy.alt), length(nid.tip))
}

test.depthTips.phylo4 <- function() {
  edgeLengthVec <- c(1.2, 1.8, 1.8, 2.1, 2.3)
  names(edgeLengthVec) <- tipLabels(phy.alt)
  checkEquals(depthTips(phy.alt), edgeLengthVec)
  tmpPhy <- phy.alt
  edgeLength(tmpPhy) <- NA
  checkTrue(is.null(depthTips(tmpPhy)))
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

test.nodeDepth.phylo4 <- function() {
  allDepths <- c(1.2, 1.8, 1.8, 2.1, 2.3, 0.9, 1.0, 1.2, 1.6)
  names(allDepths) <- names(getNode(phy.alt))
  checkIdentical(nodeDepth(phy.alt), allDepths)
  checkIdentical(nodeDepth(phy.alt, 1), allDepths[1])
  checkIdentical(nodeDepth(phy.alt, "t1"), allDepths[1])
  tmpPhy <- phy.alt
  edgeLength(tmpPhy) <- NA
  checkTrue(is.null(nodeDepth(tmpPhy)))
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
  checkIdentical(edgeLength(phy.alt), setNames(elen, eid))
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
  ## dropping all should produce empty slot
  edgeLength(phy.alt) <- numeric()
  checkIdentical(edgeLength(phy.alt), setNames(rep(NA_real_, 9), eid))
  checkIdentical(phy.alt@edge.length, numeric())
  edgeLength(phy.alt) <- NA_real_
  checkIdentical(edgeLength(phy.alt), setNames(rep(NA_real_, 9), eid))
  checkIdentical(phy.alt@edge.length, numeric())

  #
  # complete replacement
  #

  # vector with reversed names, which get matched by default
  edgeLength(phy.alt) <- numeric()
  edgeLength(phy.alt) <- setNames(elen, rev(eid))
  checkIdentical(edgeLength(phy.alt), setNames(rev(elen), eid))
  # vector with reversed names, but specify no matching
  edgeLength(phy.alt) <- numeric()
  edgeLength(phy.alt, use.names=FALSE) <- setNames(elen, rev(eid))
  checkIdentical(edgeLength(phy.alt), setNames(elen, eid))
  # vector with no names, should match to edgeId order
  edgeLength(phy.alt) <- numeric()
  edgeLength(phy.alt) <- elen
  checkIdentical(edgeLength(phy.alt), setNames(elen, eid))

  # recycling applies if fewer the nEdges elements are supplied
  # (duplicate edge length are okay)
  edgeLength(phy.alt) <- 1
  checkIdentical(edgeLength(phy.alt), setNames(rep(1, 9), eid))

  #
  # partial replacement
  #

  edgeLength(phy.alt) <- elen
  # replace an edge length using numeric index
  edgeLength(phy.alt)[9] <- 83
  checkIdentical(edgeLength(phy.alt), setNames(c(elen[1:8], 83), eid))
  # and back again, now using character index
  edgeLength(phy.alt)["8-3"] <- 0.3
  checkIdentical(edgeLength(phy.alt), setNames(elen, eid))
  # error to add length for edges that don't exist
  checkException(edgeLength(phy.alt)["fake"] <- 999)
  checkException(edgeLength(phy.alt)[999] <- 999)
  # NAs permitted only for root edge (or for *all* edges)
  edgeLength(phy.alt)[edgeId(phy.alt, "root")] <- NA
  checkIdentical(edgeLength(phy.alt), setNames(c(NA, elen[2:9]), eid))
  checkException(edgeLength(phy.alt)["8-3"] <- NA)
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

  ## dropping all should produce default tip labels, no internal labels
  labels(phy.alt) <- character()
  checkIdentical(labels(phy.alt), setNames(c(paste("T", 1:5, sep=""),
      rep(NA, 4)), nid.all))

  #
  # complete replacement
  #

  # vector with reversed names, but names not used
  labels(phy.alt) <- character()
  labels(phy.alt) <- setNames(lab.all, rev(nid.all))
  checkIdentical(labels(phy.alt), setNames(lab.all, nid.all))
  labels(phy.alt) <- character()
  labels(phy.alt, "tip") <- setNames(lab.tip, rev(nid.tip))
  checkIdentical(tipLabels(phy.alt), setNames(lab.tip, nid.tip))
  labels(phy.alt) <- character()
  labels(phy.alt, "internal") <- setNames(lab.int, rev(nid.int))
  checkIdentical(nodeLabels(phy.alt), setNames(lab.int, nid.int))
  # as above, but specify name matching, hence labels get reversed too
  labels(phy.alt) <- character()
  labels(phy.alt, use.names=TRUE) <- setNames(lab.all, rev(nid.all))
  checkIdentical(labels(phy.alt), setNames(rev(lab.all), nid.all))
  labels(phy.alt) <- character()
  labels(phy.alt, "tip", use.names=TRUE) <- setNames(lab.tip, rev(nid.tip))
  checkIdentical(tipLabels(phy.alt), setNames(rev(lab.tip), nid.tip))
  labels(phy.alt) <- character()
  labels(phy.alt, "internal", use.names=TRUE) <- setNames(lab.int, rev(nid.int))
  checkIdentical(nodeLabels(phy.alt), setNames(rev(lab.int), nid.int))
  # vector with no names, should match to nodeId order
  labels(phy.alt) <- character()
  labels(phy.alt) <- lab.all
  checkIdentical(labels(phy.alt), setNames(lab.all, nid.all))
  labels(phy.alt) <- character()
  labels(phy.alt, type="tip") <- lab.tip
  checkIdentical(tipLabels(phy.alt), setNames(lab.tip, nid.tip))
  labels(phy.alt) <- character()
  labels(phy.alt, type="internal") <- lab.int
  checkIdentical(nodeLabels(phy.alt), setNames(lab.int, nid.int))

  #
  # partial replacement
  #

  labels(phy.alt) <- lab.all
  # replace a tip using numeric index
  labels(phy.alt)[5] <- "t5a"
  checkIdentical(tipLabels(phy.alt), setNames(c(lab.tip[1:4], "t5a"), nid.tip))
  # and back again, now using character index
  labels(phy.alt)["5"] <- "t5"
  checkIdentical(labels(phy.alt), setNames(lab.all, nid.all))
  # replace an internal node using numeric index
  labels(phy.alt)[9] <- "n9a"
  checkIdentical(nodeLabels(phy.alt), setNames(c(lab.int[1:3], "n9a"), nid.int))
  # and back again, now using character index
  labels(phy.alt)["9"] <- "n9"
  checkIdentical(labels(phy.alt), setNames(lab.all, nid.all))
  # error to produce duplicate tip or internal label
  phylobase.options(allow.duplicated.labels="fail")
  checkException(labels(phy.alt)[1] <- "t2")
  checkException(labels(phy.alt)[6] <- "n7")
  # no error in allow.duplicated.labels is ok
  phylobase.options(allow.duplicated.labels="ok")
  labels(phy.alt)[1] <- "t2"
  labels(phy.alt)[6] <- "n7"
  checkIdentical(tipLabels(phy.alt), setNames(c("t2", "t2", "t3", "t4", "t5"), nid.tip))
  checkIdentical(nodeLabels(phy.alt), setNames(c("n7", "n7", "n8", "n9"), nid.int))
  # error to add labels for nodes that don't exist
  checkException(labels(phy.alt)["fake"] <- "xxx")
  checkException(labels(phy.alt)[999] <- "xxx")

}

test.nodeLabels.phylo4 <- function() {
  checkIdentical(nodeLabels(phy.alt), setNames(lab.int, nid.int))
}

test.hasNodeLabels.phylo4 <- function() {
  checkTrue(hasNodeLabels(phy.alt))
  nodeLabels(phy.alt) <- NA_character_
  checkTrue(!hasNodeLabels(phy.alt))
}

test.Replace.nodeLabels.phylo4 <- function() {

  ## dropping all should produce no internal labels
  nodeLabels(phy.alt) <- character()
  checkTrue(!any(nid.int %in% names(phy.alt@label)))
  checkIdentical(nodeLabels(phy.alt), setNames(rep(NA_character_, 4), nid.int))

  #
  # partial replacement
  #

  labels(phy.alt) <- lab.all
  # replace an internal node using numeric index
  nodeLabels(phy.alt)[4] <- "n9a"
  checkIdentical(nodeLabels(phy.alt), setNames(c(lab.int[1:3], "n9a"), nid.int))
  # and back again, now using character index
  nodeLabels(phy.alt)["9"] <- "n9"
  checkIdentical(labels(phy.alt), setNames(lab.all, nid.all))
  # error to produce duplicate internal label
  phylobase.options(allow.duplicated.labels="fail")
  checkException(nodeLabels(phy.alt)["6"] <- "n7")
  phylobase.options(op)
  phylobase.options(allow.duplicated.labels="ok")
  nodeLabels(phy.alt)["6"] <- "n7"
  checkIdentical(nodeLabels(phy.alt), setNames(c("n7", "n7", "n8", "n9"), nid.int))
  phylobase.options(op)
  # error to add labels for nodes that don't exist
  checkException(nodeLabels(phy.alt)["fake"] <- "xxx")
  checkException(nodeLabels(phy.alt)[999] <- "xxx")
}

test.tipLabels.phylo4 <- function() {
  checkIdentical(tipLabels(phy.alt), setNames(lab.tip, nid.tip))
}

test.Replace.tipLabels.phylo4 <- function() {

  ## dropping all tip labels should produce default labels
  tipLabels(phy.alt) <- character()
  checkIdentical(tipLabels(phy.alt), setNames(paste("T", 1:5, sep=""), nid.tip))

  #
  # partial replacement
  #

  labels(phy.alt) <- lab.all
  # replace a tip using numeric index
  tipLabels(phy.alt)[5] <- "t5a"
  checkIdentical(tipLabels(phy.alt), setNames(c(lab.tip[1:4], "t5a"), nid.tip))
  # and back again, now using character index
  tipLabels(phy.alt)["5"] <- "t5"
  checkIdentical(labels(phy.alt), setNames(lab.all, nid.all))
  # error to produce duplicate tip or internal label
  phylobase.options(allow.duplicated.labels="fail")
  checkException(tipLabels(phy.alt)[1] <- "t2")
  phylobase.options(op)
  phylobase.options(allow.duplicated.labels="ok")
  tipLabels(phy.alt)[1] <- "t2"
  checkIdentical(tipLabels(phy.alt), setNames(c("t2", "t2", "t3", "t4", "t5"), nid.tip))
  phylobase.options(op)
  # error to add labels for nodes that don't exist
  checkException(tipLabels(phy.alt)["fake"] <- "xxx")
  checkException(tipLabels(phy.alt)[999] <- "xxx")
}

test.hasEdgeLabels.phylo4 <- function() {
  checkTrue(hasEdgeLabels(phy.alt))
  phy.alt@edge.label <- NA_character_
  checkTrue(!hasEdgeLabels(phy.alt))
}

test.edgeLabels.phylo4 <- function() {

  # basic usage
  checkIdentical(edgeLabels(phy.alt), setNames(elab, eid))
  # should return named vector of NAs if edge labels are missing or NA
  phy.alt@edge.label <- NA_character_
  checkIdentical(edgeLabels(phy.alt), setNames(rep(NA_character_, 9), eid))
  phy.alt@edge.label <- character()
  checkIdentical(edgeLabels(phy.alt), setNames(rep(NA_character_, 9), eid))
  # if only some labels exists, should fill in NA for the others
  phy.alt@edge.label <- setNames(elab[-1], eid[-1])
  checkIdentical(edgeLabels(phy.alt), setNames(c(NA, elab[-1]), eid))

}

test.Replace.edgeLabels.phylo4 <- function() {

  ## dropping all should produce empty slot
  edgeLabels(phy.alt) <- character()
  checkIdentical(edgeLabels(phy.alt), setNames(rep(NA_character_, 9), eid))

  #
  # complete replacement
  #

  # vector with reversed names, which always get matched
  edgeLabels(phy.alt) <- character()
  edgeLabels(phy.alt) <- setNames(elab, rev(eid))
  checkIdentical(edgeLabels(phy.alt), setNames(rev(elab), eid))
  # vector with no names, should match to edgeId order
  edgeLabels(phy.alt) <- character()
  edgeLabels(phy.alt) <- elab
  checkIdentical(edgeLabels(phy.alt), setNames(elab, eid))

  # recycling applies if fewer the nEdges elements are supplied
  # (duplicate edge labels are okay)
  edgeLabels(phy.alt) <- "x"
  checkIdentical(edgeLabels(phy.alt), setNames(rep("x", 9), eid))

  #
  # partial replacement
  #

  edgeLabels(phy.alt) <- elab
  # replace an edge label using numeric index
  edgeLabels(phy.alt)[9] <- "e8-3a"
  checkIdentical(edgeLabels(phy.alt), setNames(c(elab[1:8], "e8-3a"), eid))
  # and back again, now using character index
  edgeLabels(phy.alt)["8-3"] <- "e8-3"
  checkIdentical(edgeLabels(phy.alt), setNames(elab, eid))
  # error to add labels for edges that don't exist
  checkException(edgeLabels(phy.alt)["fake"] <- "xxx")
  checkException(edgeLabels(phy.alt)[999] <- "xxx")
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

test.isUltrametric <- function() {
  checkTrue(!isUltrametric(phy.alt))
  tmpPhy <- as(rcoal(10), "phylo4")
  checkTrue(isUltrametric(tmpPhy))
  tmpPhy <- phy.alt
  edgeLength(tmpPhy) <- NA
  checkException(isUltrametric(tmpPhy))
}

phylobase.options(op)
