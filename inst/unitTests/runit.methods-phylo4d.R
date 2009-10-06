#
# --- Test methods-phylo4d.R ---
#

# create phylo4 object with a full complement of valid slots
ancestor <- as.integer(c(6,7,7,6,8,0,8,9,9))
descendant <- as.integer(c(7,1,2,8,3,6,9,4,5))
edge <- cbind(ancestor, descendant)
nid.tip <- 1:5
nid.int <- 6:9
nid.all <- c(nid.tip, nid.int)
lab.tip <- paste("t", nid.tip, sep="")
lab.int <- paste("n", nid.int, sep="")
lab.all <- c(lab.tip, lab.int)
elen <- descendant/10
elab <- paste("e", ancestor, descendant, sep="-")
phy <- phylo4(x=edge, tip.label=lab.tip, node.label=lab.int,
    edge.length=elen, edge.label=elab)

# now create phylo4d by adding data (with node IDs as row.names)
allDt <- data.frame(a=letters[nid.all], b=10*nid.all)
tipDt <- data.frame(c=letters[nid.tip], d=10*nid.tip)
nodDt <- data.frame(c=letters[nid.int], e=10*nid.int)
row.names(allDt) <- nid.all
row.names(tipDt) <- nid.tip
row.names(nodDt) <- nid.int
phyd <- phylo4d(phy, tip.data=tipDt, node.data=nodDt, all.data=allDt,
                match.data=TRUE, merge.data=TRUE)

# create altered version such that each slot is out of order with
# respect to all others; methods should be able to handle this
phyd.alt <- phyd
phyd.alt@label <- rev(phyd@label)
phyd.alt@edge <- phyd@edge[c(6:9, 1:5), ]
phyd.alt@edge.length <- phyd@edge.length[c(7:9, 1:6)]
phyd.alt@edge.label <- phyd@edge.label[c(8:9, 1:7)]
nid.tip.r <- c(2,5,4,3,1)
nid.int.r <- c(8,7,9,6)
nid.all.r <- c(nid.tip.r, nid.int.r)
phyd.alt@data <- phyd@data[rank(nid.all.r), ]

# for comparisons, manually create expected "all" trait data.frame
m1 <- merge(allDt, rbind(tipDt["c"], nodDt["c"]), by=0, all=TRUE)
m2 <- merge(tipDt["d"], nodDt["e"], by=0, all=TRUE)
eAllDt <- merge(m1, m2, by="Row.names", all=TRUE)[-1]
row.names(eAllDt) <- lab.all

# for comparisons, manually create expected "tip" trait data.frame
m1 <- merge(allDt, rbind(tipDt["c"], nodDt["c"]), by=0, all=TRUE)
m2 <- merge(tipDt["d"], nodDt["e"], by=0, all=TRUE)
eTipDt <- merge(m1, m2, by="Row.names", all=TRUE)[nid.tip, -1]
row.names(eTipDt) <- lab.tip

# manually create expected tip trait data.frame
m1 <- merge(allDt, rbind(tipDt["c"], nodDt["c"]), by=0, all=TRUE)
m2 <- merge(tipDt["d"], nodDt["e"], by=0, all=TRUE)
eNodDt <- merge(m1, m2, by="Row.names", all=TRUE)[nid.int, -1]
row.names(eNodDt) <- lab.int

#-----------------------------------------------------------------------

test.tdata.phylo4d <- function() {
    # function(x, type=c("tip", "internal", "allnode"),
    #   label.type=c("row.names","column"), empty.columns=TRUE, ...)

    # check basic tdata usage
    checkIdentical(tdata(phyd.alt, type="tip"), eTipDt)
    checkIdentical(tdata(phyd.alt, type="internal"), eNodDt)
    checkIdentical(tdata(phyd.alt, type="all"), eAllDt)

    # label.type="row.names"
    tmpDt <- data.frame(eAllDt[nid.tip, -5, ], row.names=lab.tip)
    checkIdentical(tdata(phyd.alt, type="tip", label.type="row.names",
        empty.columns=FALSE), data.frame(tmpDt[nid.tip,], row.names=lab.tip))
    # label.type="column"
    tmpDt <- data.frame(label=lab.tip, eAllDt[nid.tip, -5, ],
        row.names=as.character(nid.tip))
    checkIdentical(tdata(phyd.alt, type="tip", label.type="column",
                         empty.columns=FALSE), tmpDt)

    # keep empty.columns
    checkIdentical(tdata(phyd.alt, type="tip", empty.columns=TRUE),
        eAllDt[nid.tip,])

    #
    # misc tests
    #

    # check with other tree orderings
    phyd.pre <- reorder(phyd.alt, "preorder")
    checkIdentical(tdata(phyd.pre, "all", empty.columns=FALSE), eAllDt)
    phyd.post <- reorder(phyd.alt, "postorder")
    checkIdentical(tdata(phyd.post, "all", empty.columns=FALSE), eAllDt)

}

## currently just basic tests of tdata replacement; using out-of-order
## data, but only with default args (e.g. row.name-nodeID matching)
## ... formatData unit tests should be sufficient for the rest
test.Replace.tdata.phylo4d <- function() {

    ## replace data, labels are row names
    tdata(phyd.alt, "all") <- allDt[rank(nid.all.r), , drop=FALSE]
    checkIdentical(tdata(phyd.alt, type="all"), data.frame(allDt,
        row.names=lab.all))

    ## replace data with empty data frame
    tdata(phyd.alt) <- data.frame()
    checkIdentical(tdata(phyd.alt), data.frame(row.names=lab.all))

    ## same as first test, but leaving out default 'all' type
    tdata(phyd.alt) <- allDt[rank(nid.all.r), , drop=FALSE]
    checkIdentical(tdata(phyd.alt), data.frame(allDt,
        row.names=lab.all))

}

test.tipData.phylo4d <- function() {
    # label.type="row.names"
    checkIdentical(tipData(phyd.alt, label.type="row.names",
        empty.columns=FALSE), eTipDt[-5])
    # label.type="column"
    tmpDt <- data.frame(label=lab.tip, eTipDt[-5],
        row.names=as.character(nid.tip))
    checkIdentical(tipData(phyd.alt, label.type="column",
        empty.columns=FALSE), tmpDt)

    # keep empty.columns
    checkIdentical(tipData(phyd.alt), eTipDt)
}

test.Replace.tipData.phylo4d <- function() {
    ## replace data with tip data only, clearing all data
    tipData(phyd.alt, clear.all=TRUE) <- tipDt[rank(nid.tip.r), ,
        drop=FALSE]
    checkIdentical(tipData(phyd.alt), data.frame(tipDt,
        row.names=lab.tip))
}

test.nodeData.phylo4d <- function() {

    # label.type="row.names"
    checkIdentical(nodeData(phyd.alt, label.type="row.names",
        empty.columns=FALSE), eNodDt[-4])

    # label.type="column"
    tmpDt <- data.frame(label=lab.int, eNodDt[-4],
        row.names=as.character(nid.int))
    checkIdentical(nodeData(phyd.alt, label.type="column",
        empty.columns=FALSE), tmpDt)

    # keep empty.columns
    checkIdentical(nodeData(phyd.alt), eNodDt)
}

test.Replace.nodeData.phylo4d <- function() {
    ## replace data with internal data only, clearing all data
    nodeData(phyd.alt, clear.all=TRUE) <- nodDt[rank(nid.int.r), ,
        drop=FALSE]
    checkIdentical(nodeData(phyd.alt), data.frame(nodDt,
        row.names=lab.int))
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

