#
# --- Test subset.R ---
#

# load test comparison objects
load("trees.RData")
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((t1:0.2,(t2:0.1,t3:0.1):0.15):0.5,t4:0.7):0.2,t5:1):0.4;")
tr.sub2 <- read.tree(text="(t2:0.95,t5:1);")
tr.sub4 <- read.tree(text="(((t1:0.2,t2:0.25):0.5,t4:0.7):0.2,t5:1);")
# Explicitly set order as 'cladewise', to match behavior of
# ape::drop.tip in test comparisons below
tr <- reorder(tr, "cladewise")
tr.sub2 <- reorder(tr.sub2, "cladewise")
tr.sub4 <- reorder(tr.sub4, "cladewise")

test.subset.phylo <- function() {
    # subset 2 tips
    checkEquals(tr.sub2, subset(tr, tips.include=c(2, 5)))
    checkEquals(tr.sub2, subset(tr, tips.exclude=c(1, 3, 4)))
    checkEquals(tr.sub2, subset(tr, tips.include=c("t2", "t5")))
    checkEquals(tr.sub2, subset(tr, tips.exclude=c("t1", "t3", "t4")))
    # subset 4 tips
    checkEquals(tr.sub4, subset(tr, tips.include=c(1, 2, 4, 5)))
    checkEquals(tr.sub4, subset(tr, tips.exclude=3))
    checkEquals(tr.sub4, subset(tr, tips.include=c("t1", "t2", "t4", "t5")))
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
    checkException(subset(phy, tips.include="t1"))
    checkException(subset(phy, tips.include=c("t1", "t999")))
    # error if zero valid tips requested
    checkException(subset(phy, tips.include="t999"))
    # error if more than one subset criteria are supplied
    checkException(subset(phyd, tips.include="t1", tips.exclude="t3"))
}

test.subset.phylo4d <- function() {
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
    # extract phylo4 from itself
    phy <- phylo4(tr, annote=list(x="annotation"))
    checkIdentical(phy, extractTree(phy))

    # extract phylo4 from phylo4d
    phyd <- phylo4d(tr, tip.data= data.frame(x=1:5, row.names=tr$tip.label),
      annote=list(x="annotation"), metadata=list(x="metadata"))
    checkIdentical(phy, extractTree(phyd))
}
