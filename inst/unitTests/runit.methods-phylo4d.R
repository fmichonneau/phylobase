#
# --- Test methods-phylo4.R ---
#

# Create sample tree for testing (ape::phylo object)
tr <- read.tree(text="(((spA:0.2,(spB:0.1,spC:0.1):0.15):0.5,spD:0.7):0.2,spE:1):0.4;")
phy <- as(tr, "phylo4")
##   label node ancestr edge.length node.type
## 6  <NA>    6      NA        0.40      root
## 7  <NA>    7       6        0.20  internal
## 8  <NA>    8       7        0.50  internal
## 9  <NA>    9       8        0.15  internal
## 1   spA    1       8        0.20       tip
## 2   spB    2       9        0.10       tip
## 3   spC    3       9        0.10       tip
## 4   spD    4       7        0.70       tip
## 5   spE    5       6        1.00       tip
allDt <- data.frame(a=rnorm(nTips(phy)+nNodes(phy)))
tipDt <- data.frame(b=letters[1:nTips(phy)], c=rnorm(nTips(phy)))
nodDt <- data.frame(d=rnorm(nNodes(phy)))
rownames(tipDt) <- 1:nTips(phy)
rownames(nodDt) <- (nTips(phy)+1):(nTips(phy)+nNodes(phy))
rownames(allDt) <- 1:(nTips(phy)+nNodes(phy))
phyd <- phylo4d(phy, tip.data=tipDt, node.data=nodDt, all.data=allDt,
                match.data=TRUE, merge.data=TRUE)


test.tdata.phylo4d <- function() {
    # function(x, type=c("tip", "internal", "allnode"),
    #   label.type=c("row.names","column"), empty.columns=TRUE, ...)
    compDt <- cbind(tipLabels(phy), a=allDt[1:nTips(phy), ], tipDt)
    colnames(compDt)[1] <- "label"
    checkIdentical(tdata(phyd, type="tip", label.type="column",
                         empty.columns=FALSE), compDt)
    phyd <- reorder(phyd, "preorder")
    checkIdentical(tdata(phyd, label.type="column",
                         empty.columns=FALSE), compDt)
    phyd <- reorder(phyd, "postorder")
    checkIdentical(tdata(phyd, label.type="column",
                         empty.columns=FALSE), compDt)

}

test.Replace.tdata.phylo4d <- function() {
    # function(object, type = c("tip", "internal", "allnode"), ...,
    #   value)

    op <- options()
    options(stringsAsFactors=FALSE)

    newTipDt <- data.frame(a=1:nTips(phy), b=10+(1:nTips(phy)))
    rownames(newTipDt) <- nodeId(phy, "tip")
    newNodDt <- data.frame(c=1:nNodes(phy), d=10+(1:nNodes(phy)))
    rownames(newNodDt) <- nodeId(phy, "internal")
    newAllDt <- data.frame(e=1:(nTips(phy)+nNodes(phy)))
    rownames(newAllDt) <- nodeId(phy, "all")

    ## default ordering
    tmpTip <- tmpNod <- tmpAll <- phyd
    tdata(tmpTip, type="tip") <- newTipDt
    compTip <- data.frame(label=tipLabels(tmpTip), newTipDt)
    checkIdentical(tdata(tmpTip, type="tip", label.type="column"), compTip)

    tdata(tmpNod, type="internal") <- newNodDt
    compNod <- data.frame(label=nodeId(tmpNod, "internal"), newNodDt)
    rownames(compNod) <- as.character(nodeId(tmpNod, "internal"))
    checkIdentical(tdata(tmpNod, type="internal", label.type="column",
                         empty.columns=FALSE), compNod)

    tdata(tmpAll, type="all") <- newAllDt
    compAll <- data.frame(label=c(tipLabels(tmpNod), nodeId(tmpNod, "internal")),
                          newAllDt)
    checkIdentical(tdata(tmpAll, type="all", label.type="column"), compAll)

    ## Preorder
    tmpTip <- tmpNod <- tmpAll <- reorder(phyd, "preorder")
    tdata(tmpTip, type="tip") <- newTipDt
    compTip <- data.frame(label=tipLabels(tmpTip), newTipDt)
    checkIdentical(tdata(tmpTip, type="tip", label.type="column"), compTip)

    tdata(tmpNod, type="internal") <- newNodDt
    compNod <- data.frame(label=nodeId(tmpNod, "internal"), newNodDt)
    rownames(compNod) <- as.character(nodeId(tmpNod, "internal"))
    checkIdentical(tdata(tmpNod, type="internal", label.type="column",
                         empty.columns=FALSE), compNod)

    tdata(tmpAll, type="all") <- newAllDt
    compAll <- data.frame(label=c(tipLabels(tmpNod), nodeId(tmpNod, "internal")),
                          newAllDt)
    checkIdentical(tdata(tmpAll, type="all", label.type="column"), compAll)

    ## Postorder
    tmpTip <- tmpNod <- tmpAll <- reorder(phyd, "postorder")
    tdata(tmpTip, type="tip") <- newTipDt
    compTip <- data.frame(label=tipLabels(tmpTip), newTipDt)
    checkIdentical(tdata(tmpTip, type="tip", label.type="column"), compTip)

    tdata(tmpNod, type="internal") <- newNodDt
    compNod <- data.frame(label=nodeId(tmpNod, "internal"), newNodDt)
    rownames(compNod) <- as.character(nodeId(tmpNod, "internal"))
    checkIdentical(tdata(tmpNod, type="internal", label.type="column"), compNod)

    tdata(tmpAll, type="all") <- newAllDt
    compAll <- data.frame(label=c(tipLabels(tmpNod), nodeId(tmpNod, "internal")),
                          newAllDt)
    checkIdentical(tdata(tmpAll, type="all", label.type="column"), compAll)

    options(op)
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

