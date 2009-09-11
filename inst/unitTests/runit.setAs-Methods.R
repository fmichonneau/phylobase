#
# --- Test setAs-Methods.R ---
#

# create ape::phylo version of a simple tree for testing
nwk <- "((t1:0.1,t2:0.2)n7:0.7,(t3:0.3,(t4:0.4,t5:0.5)n9:0.9)n8:0.8)n6:0.6;"
tr <- read.tree(text=nwk)

# create analogous phylo4 object with a full complement of valid slots
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

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phy.alt <- phy
phy.alt@tip.label <- rev(phy@tip.label)
phy.alt@node.label <- rev(phy@node.label)
phy.alt@edge <- phy@edge[c(6:9, 1:5), ]
phy.alt@edge.length <- phy@edge.length[c(7:9, 1:6)]
phy.alt@edge.label <- phy@edge.label[c(8:9, 1:7)]

#-----------------------------------------------------------------------
 
test.phylo.As.phylo4 <- function() {
  checkIdentical(as(tr, "phylo4"), phylo4(tr))
}

test.phylo.As.phylo4d <- function() {
  checkIdentical(as(tr, "phylo4d"), phylo4d(tr))
}

test.multiPhylo.As.multiPhylo4 <- function() {
}

test.multiPhylo4.As.multiPhylo <- function() {
}

test.phylo4.As.phylo <- function() {
# note: checkEquals("phylo") uses all.equal.phylo()

  # phylo tree in unknown order
  checkEquals(suppressWarnings(as(phy, "phylo")), tr)
  # ...now check for warning for unknown order
  opt <- options(warn=3)
  checkException(as(phy, "phylo"))
  options(opt)

  # phylo tree in cladewise order
  tr.cladewise <- reorder(tr, "cladewise")
  phy.c <- as(tr.cladewise, "phylo4")
  checkEquals(as(phy.c, "phylo"), tr.cladewise)

  # phylo tree in pruningwise order
  tr.pruningwise <- reorder(tr, "pruningwise")
  phy.p <- as(tr.pruningwise, "phylo4")
  checkEquals(as(phy.p, "phylo"), tr.pruningwise)

  # after transforming the jumbled tree to phylo and back, edge matrix
  # and edge slots should still be in the original order, but node slots
  # should be back in nodeId order
  phy.r <- reorder(phy.alt)
  phy.roundtrip.r <- reorder(phylo4(as(phy.alt, "phylo")))
  checkIdentical(edges(phy.roundtrip.r), edges(phy.r))
  checkIdentical(edgeLength(phy.roundtrip.r), edgeLength(phy.r))
  checkIdentical(labels(phy.roundtrip.r), labels(phy.r))

}

# this coerce method is defined implicitly
test.phylo4d.As.phylo <- function() {
# note: checkEquals("phylo") uses all.equal.phylo()

  # phylo tree in unknown order
  phyd <- as(tr, "phylo4d")
  tdata(phyd) <- data.frame(x=1:5, row.names=tipLabels(phyd))
  checkEquals(as(phyd, "phylo"), tr)
  # ...now check for warning for unknown order
  opt <- options(warn=3)
  checkException(as(phyd, "phylo"))
  options(opt)

  # phylo tree in cladewise order
  tr.cladewise <- reorder(tr, "cladewise")
  phyd <- as(tr.cladewise, "phylo4d")
  tdata(phyd) <- data.frame(x=1:5, row.names=tipLabels(phyd))
  checkEquals(as(phyd, "phylo"), tr.cladewise)
  # ...now check for warning for dropping data
  opt <- options(warn=3)
  checkException(as(phyd, "phylo"))
  options(opt)

  # phylo tree in pruningwise order
  tr.pruningwise <- reorder(tr, "pruningwise")
  phyd <- as(tr.pruningwise, "phylo4d")
  tdata(phyd) <- data.frame(x=1:5, row.names=tipLabels(phyd))
  checkEquals(as(phyd, "phylo"), tr.pruningwise)
}

test.phylo4.As.phylog <- function() {
}

test..phylo4ToDataFrame <- function() {
  phy.show <- phylobase:::.phylo4ToDataFrame(phy.alt, "pretty")
  checkIdentical(phy.show$label, c(lab.tip, lab.int))
  checkIdentical(phy.show$node, c(nid.tip, nid.int))
  checkIdentical(phy.show$ancestor, ancestor[match(c(nid.tip, nid.int),
    descendant)])
  checkIdentical(phy.show$edge.length, sort(elen))
  checkIdentical(phy.show$node.typ, factor(unname(nodeType(phy))))
}

## core functionality is already tested in test..phylo4ToDataFrame()
test.phylo4.As.data.frame <- function() {
    # rooted tree
    checkTrue(is.data.frame(as(phy, "data.frame")))

    # unrooted tree
    tru <- unroot(tr)
    phyu <- as(tru, "phylo4")
    # should probably check that this coercion results in something
    # *correct*, not just that it produces a data.frame
    checkTrue(is.data.frame(as(phyu, "data.frame")))
}
