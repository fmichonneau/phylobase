#
# --- Test subset.R ---
#
 
# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;") 

test.subset.phylo <- function() {
    print(subset(tr, 1:4))
    print(subset(tr, 1:2))
    checkEquals(tr, subset(tr))
}

test.subset.phylo4 <- function() {
    phy <- as(tr, "phylo4")
    print(subset(phy, 1:4))
    print(subset(phy, 1:2))
    # check variants that should all return the original object
    checkEquals(phy, subset(phy))
    checkEquals(phy, subset(phy, tipLabels(phy)))
    checkEquals(phy, subset(phy, seq_len(nTips(phy))))
    checkEquals(phy, phy[tipLabels(phy)])
    checkEquals(phy, phy[seq_len(nTips(phy))])
    checkException(subset(phy,"spA"))
    checkException(subset(phy,c("spA","spZ")))
    checkException(subset(phy,c("spZ")))
}

test.subset.phylo4d <- function() {
    phyd <- as(tr, "phylo4d")
    tdata(phyd) <- 1:5
    print(subset(phyd, 1:4))
    print(subset(phyd, 1:2))
    # check variants that should all return the original object
    checkEquals(phyd, subset(phyd))
    checkEquals(phyd, subset(phyd, tipLabels(phyd)))
    checkEquals(phyd, subset(phyd, seq_len(nTips(phyd))))
    checkEquals(phyd, phyd[tipLabels(phyd)])
    checkEquals(phyd, phyd[seq_len(nTips(phyd))])
    set.seed(1)
    z1 = rcoal(5)
    zdat = data.frame(1:4)
    rownames(zdat) <- paste("t",1:4,sep="")
    z2 = phylo4d(z1,zdat,missing="OK")
    subset(z2,tips.exclude="t1")
    subset(z2,tips.include=c("t4","t5"))
}

test.extractTree <- function() {
    phyd <- as(tr, "phylo4d")
    phy <- as(tr, "phylo4")
    checkEquals(phy, extractTree(phyd))
}
