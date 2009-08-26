#
# --- Test setAs-Methods.R ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 
phy <- as(tr, "phylo4")

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
  phy <- as(tr, "phylo4")
  checkEquals(as(phy, "phylo"), tr)
  # ...now check for warning for unknown order
  opt <- options(warn=3)
  checkException(as(phy, "phylo"))
  options(opt)

  # phylo tree in cladewise order
  tr.cladewise <- reorder(tr, "cladewise")
  phy <- as(tr.cladewise, "phylo4")
  checkEquals(as(phy, "phylo"), tr.cladewise)

  # phylo tree in pruningwise order
  tr.pruningwise <- reorder(tr, "pruningwise")
  phy <- as(tr.pruningwise, "phylo4")
  checkEquals(as(phy, "phylo"), tr.pruningwise)
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
