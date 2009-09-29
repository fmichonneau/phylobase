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
phyd.alt@tip.label <- rev(phyd@tip.label)
phyd.alt@node.label <- rev(phyd@node.label)
phyd.alt@edge <- phyd@edge[c(6:9, 1:5), ]
phyd.alt@edge.length <- phyd@edge.length[c(7:9, 1:6)]
phyd.alt@edge.label <- phyd@edge.label[c(8:9, 1:7)]
nid.tip.r <- c(2,5,4,3,1)
nid.int.r <- c(8,7,9,6)
nid.all.r <- c(nid.tip.r, nid.int.r)
phyd.alt@tip.data <- phyd@tip.data[rank(nid.tip.r), ]
phyd.alt@node.data <- phyd@node.data[rank(nid.int.r), ]

#-----------------------------------------------------------------------

test.tdata.phylo4d <- function() {
DEACTIVATED("turned off until tdata can handle out-of-order rows in data slots")
    # function(x, type=c("tip", "internal", "allnode"),
    #   label.type=c("row.names","column"), empty.columns=TRUE, ...)

    # TODO: flesh out these tests!
    tip.data <- tdata(phyd.alt, type="tip")
    checkTrue(is.data.frame(tip.data))
    int.data <- tdata(phyd.alt, type="internal")
    checkTrue(is.data.frame(int.data))
    all.data <- tdata(phyd.alt, type="all")
    checkTrue(is.data.frame(all.data))

    #
    # label.type
    #

    # label.type="row.names"
    compDt <- data.frame(allDt[nid.tip, ], tipDt, row.names=lab.tip)
    checkIdentical(tdata(phyd.alt, type="tip", label.type="row.names",
                         empty.columns=FALSE), compDt)
    # label.type="column"
    compDt <- cbind(label=lab.tip, allDt[nid.tip, ], tipDt,
        row.names=as.character(nid.tip))
    checkIdentical(tdata(phyd.alt, type="tip", label.type="column",
                         empty.columns=FALSE), compDt)

    #
    # empty.columns
    #

    #
    # misc tests
    #

    # check with other tree orderings
    phyd.pre <- reorder(phyd.alt, "preorder")
    checkIdentical(tdata(phyd.pre, label.type="column",
                         empty.columns=FALSE), compDt)
    phyd.post <- reorder(phyd.alt, "postorder")
    checkIdentical(tdata(phyd.post, label.type="column",
                         empty.columns=FALSE), compDt)

}

## currently just basic tests of tdata replacement; using out-of-order
## data, but only with default args (e.g. row.name-nodeID matching)
test.Replace.tdata.phylo4d <- function() {

    ## replace data with tip data only
    tdata(phyd.alt, type="tip") <- tipDt[rank(nid.tip.r), , drop=FALSE]
    checkIdentical(tdata(phyd.alt, type="tip"), data.frame(tipDt,
        row.names=lab.tip))

    ## replace data with internal data only
    tdata(phyd.alt, type="internal") <- nodDt[rank(nid.int.r), , drop=FALSE]
    checkIdentical(tdata(phyd.alt, type="internal"), data.frame(nodDt,
        row.names=lab.int))

    ## replace data with both tip and internal data
    tdata(phyd.alt, type="all") <- allDt[rank(nid.all.r), , drop=FALSE]
    checkIdentical(tdata(phyd.alt, type="all"), data.frame(allDt,
        row.names=lab.all))
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

