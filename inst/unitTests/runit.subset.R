#
# --- Test subset.R ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((t1:0.2,(t2:0.1,t3:0.1):0.15):0.5,t4:0.7):0.2,t5:1):0.4;")
tr.sub2 <- read.tree(text="(t2:0.95,t5:1);")
tr.sub4 <- read.tree(text="(((t1:0.2,t2:0.25):0.5,t4:0.7):0.2,t5:1);")

test.subset.phylo <- function() {
    # subset 2 tips
    checkEquals(tr.sub2, subset(tr, tips.include=c(2, 5)))
    checkEquals(tr.sub2, subset(tr, tips.exclude=c(1, 3, 4)))
    checkEquals(tr.sub2, subset(tr, tips.include=c("t2", "t5")))
    checkEquals(tr.sub2, subset(tr, tips.exclude=c("t1", "t3", "t4")))
    # subset 4 tips
    checkEquals(tr.sub4, subset(tr, tips.include=c(1, 2, 4, 5)))
    checkEquals(tr.sub4, subset(tr, tips.exclude=3))
    checkEquals(tr.sub4, subset(tr, tips.include=c("t1", "t2", "t4",
      "t5")))
    checkEquals(tr.sub4, subset(tr, tips.exclude="t3"))
    # check variants that should all return the original object
    checkEquals(tr, subset(tr))
    # error if only one valid tip requested
    checkException(subset(tr, tips.include="t1"))
    checkException(subset(tr, tips.include=c("t1", "t999")))
    # error if zero valid tips requested
    checkException(subset(tr, tips.include="t999"))
}

test.subset.phylo4 <- function() {
  DEACTIVATED("Broken?: subset changes phy order from 'unknown' to 'preorder'")
    phy <- phylo4(tr)
    phy.sub2 <- phylo4(tr.sub2)
    phy.sub4 <- phylo4(tr.sub4)
    # subset 2 tips
    checkEquals(phy.sub2, subset(phy, tips.include=c(2, 5)))
    checkEquals(phy.sub2, subset(phy, tips.exclude=c(1, 3, 4)))
    checkEquals(phy.sub2, subset(phy, tips.include=c("t2", "t5")))
    checkEquals(phy.sub2, subset(phy, tips.exclude=c("t1", "t3", "t4")))
    # subset 4 tips
    checkEquals(phy.sub4, subset(phy, tips.include=c(1, 2, 4, 5)))
    checkEquals(phy.sub4, subset(phy, tips.exclude=3))
    checkEquals(phy.sub4, subset(phy, tips.include=c("t1", "t2", "t4", "t5")))
    checkEquals(phy.sub4, subset(phy, tips.exclude="t3"))
    # check variants that should all return the original object
    checkEquals(phy, subset(phy))
    checkEquals(phy, subset(phy, tipLabels(phy)))
    checkEquals(phy, subset(phy, seq_len(nTips(phy))))
    checkEquals(phy, phy[tipLabels(phy)])
    checkEquals(phy, phy[seq_len(nTips(phy))])
    # error if only one valid tip requested
    checkException(subset(phy, tip.include="t1"))
    checkException(subset(phy, tip.include=c("t1", "t999")))
    # error if zero valid tips requested
    checkException(subset(phy, tip.include="t999"))
    # error if more than one subset criteria are supplied
    checkException(subset(phyd, tips.include="t1", tips.exclude="t3"))
}

test.subset.phylo4d <- function() {
    phyd <- phylo4d(tr, data.frame(x=1:5, row.names=paste("t", 1:5, sep="")))
    phyd.sub2 <- phylo4d(tr.sub2, data.frame(x=c(2,5),
      row.names=paste("t", c(2,5), sep="")))
    phyd.sub4 <- phylo4d(tr.sub4, data.frame(x=c(1,2,4,5),
      row.names=paste("t", c(1,2,4,5), sep="")))
    # subset 2 tips
    checkEquals(phyd.sub2, subset(phyd, tips.include=c(2, 5)))
    checkEquals(phyd.sub2, subset(phyd, tips.exclude=c(1, 3, 4)))
    checkEquals(phyd.sub2, subset(phyd, tips.include=c("t2", "t5")))
    checkEquals(phyd.sub2, subset(phyd, tips.exclude=c("t1", "t3", "t4")))
    # subset 4 tips
    checkEquals(phyd.sub4, subset(phyd, tips.include=c(1, 2, 4, 5)))
    checkEquals(phyd.sub4, subset(phyd, tips.exclude=3))
    checkEquals(phyd.sub4, subset(phyd, tips.include=c("t1", "t2", "t4", "t5")))
    checkEquals(phyd.sub4, subset(phyd, tips.exclude="t3"))
    # check variants that should all return the original object
    checkEquals(phyd, subset(phyd))
    checkEquals(phyd, subset(phyd, tipLabels(phyd)))
    checkEquals(phyd, subset(phyd, seq_len(nTips(phyd))))
    checkEquals(phyd, phyd[tipLabels(phyd)])
    checkEquals(phyd, phyd[seq_len(nTips(phyd))])
    # error if only one valid tip requested
    checkException(subset(phyd, tips.include="t1"))
    checkException(subset(phyd, tips.include=c("t1", "t999")))
    # error if zero valid tips requested
    checkException(subset(phyd, tips.include="t999"))
    # subset tips that include an NA value
    tdata(phyd)["t5", "x"] <- NA
    tdata(phyd.sub2)["t5", "x"] <- NA
    checkEquals(phyd.sub2, subset(phyd, tips.include=c(2, 5)))
    checkEquals(phyd.sub2, subset(phyd, tips.exclude=c(1, 3, 4)))
    checkEquals(phyd.sub2, subset(phyd, tips.include=c("t2", "t5")))
    checkEquals(phyd.sub2, subset(phyd, tips.exclude=c("t1", "t3", "t4")))
}

test.extractTree <- function() {
    phyd <- as(tr, "phylo4d")
    phy <- as(tr, "phylo4")
    checkEquals(phy, extractTree(phyd))
}
