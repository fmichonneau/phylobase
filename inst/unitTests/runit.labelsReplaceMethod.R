data(geospiza)

p4 <- extractTree(geospiza)
p4d <- geospiza

test.labelsTipsPhylo4 <- function() {
    tLbl <- paste("t", 1:nTips(p4), sep="")
    nmTLbl <- tLbl
    names(nmTLbl) <- nodeId(p4, "tip")
    tLbl <- sample(tLbl)
    nmTLbl <- sample(nmTLbl)

    ## case all options by default and unnamed vector
    p4c <- p4
    labels(p4c) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4c@tip.label), tLbl)

    ## case all options by default and named vector
    p4c <- p4
    labels(p4c) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4c@tip.label), unname(nmTLbl))

    ## case type defined
    p4c <- p4
    labels(p4c, "tip") <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4c@tip.label), tLbl)

    ## case type defined and use.names=TRUE but no names
    p4c <- p4
    labels(p4c, "tip", use.names=TRUE) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4c@tip.label), tLbl)

    ## case type defined and use.names=TRUE with names
    p4c <- p4
    labels(p4c, "tip", use.names=TRUE) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(p4c@tip.label, nmTLbl[order(as.numeric(names(nmTLbl)))])
}

test.labelsNodePhylo4 <- function() {

    ndLbl <- paste("n", 1:nNodes(p4), sep="")
    nmNdLbl <- ndLbl
    names(nmNdLbl) <- nodeId(p4, "internal")

    ndLbl <- sample(ndLbl)
    nmNdLbl <- sample(nmNdLbl)

    ## case type defined
    p4c <- p4
    labels(p4c, "internal") <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@node.label) %in% nodeId(p4c, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4c@node.label), ndLbl)

    ## case type defined and use.names=TRUE but no names
    p4c <- p4
    labels(p4c, "internal", use.names=TRUE) <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@node.label) %in% nodeId(p4c, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4c@node.label), ndLbl)

    ## case type defined and use.names=TRUE with names
    p4c <- p4
    labels(p4c, "internal", use.names=TRUE) <- nmNdLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4c@node.label) %in% nodeId(p4c, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(p4c@node.label, nmNdLbl[order(as.numeric(names(nmNdLbl)))])

}

test.labelsTipsPhylo4d <- function() {
    tLbl <- paste("t", 1:nTips(p4d), sep="")
    nmTLbl <- tLbl
    names(nmTLbl) <- nodeId(p4d, "tip")
    tLbl <- sample(tLbl)
    nmTLbl <- sample(nmTLbl)

    browser()
    ## case all options by default and unnamed vector
    p4dc <- p4d
    labels(p4dc) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4dc@tip.label), tLbl)

    ## case all options by default and named vector
    p4dc <- p4d
    labels(p4dc) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4dc@tip.label), unname(nmTLbl))

    ## case type defined
    p4dc <- p4d
    labels(p4dc, "tip") <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4dc@tip.label), tLbl)

    ## case type defined and use.names=TRUE but no names
    p4dc <- p4d
    labels(p4dc, "tip", use.names=TRUE) <- tLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4dc@tip.label), tLbl)

    ## case type defined and use.names=TRUE with names
    p4dc <- p4d
    labels(p4dc, "tip", use.names=TRUE) <- nmTLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(p4dc@tip.label, nmTLbl[order(as.numeric(names(nmTLbl)))])
}

test.labelsNodePhylo4d <- function() {

    ndLbl <- paste("n", 1:nNodes(p4d), sep="")
    nmNdLbl <- ndLbl
    names(nmNdLbl) <- nodeId(p4d, "internal")

    ndLbl <- sample(ndLbl)
    nmNdLbl <- sample(nmNdLbl)

    ## case type defined
    p4dc <- p4d
    labels(p4dc, "internal") <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@node.label) %in% nodeId(p4dc, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4dc@node.label), ndLbl)

    ## case type defined and use.names=TRUE but no names
    p4dc <- p4d
    labels(p4dc, "internal", use.names=TRUE) <- ndLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@node.label) %in% nodeId(p4dc, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(unname(p4dc@node.label), ndLbl)

    ## case type defined and use.names=TRUE with names
    p4dc <- p4d
    labels(p4dc, "internal", use.names=TRUE) <- nmNdLbl
    ## check the internal names are there and match tips
    checkTrue(all(names(p4dc@node.label) %in% nodeId(p4dc, "internal")))
    ## check that the labels are correct: here use.names=FALSE
    ## so the order should be the same as in the shuffled vector
    ## of name labels
    checkEquals(p4dc@node.label, nmNdLbl[order(as.numeric(names(nmNdLbl)))])

}

test.labelsAllPhylo4 <- function() {

    allLbl <- paste("n", 1:(nTips(p4)+nNodes(p4)), sep="")
    nmAllLbl <- allLbl
    names(nmAllLbl) <- nodeId(p4, "all")

    allLbl <- sample(allLbl)
    nmAllLbl <- sample(nmAllLbl)

    p4c <- p4
    labels(p4c, "all") <- allLbl
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    checkTrue(all(names(p4c@node.label) %in% nodeId(p4c, "internal")))
    checkEquals(unname(p4c@tip.label), allLbl[1:nTips(p4)])
    checkEquals(unname(p4c@node.label),
                allLbl[(nTips(p4)+1):(nTips(p4)+nNodes(p4))])

    p4c <- p4
    labels(p4c, "all") <- nmAllLbl
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    checkTrue(all(names(p4c@node.label) %in% nodeId(p4c, "internal")))
    checkEquals(unname(p4c@tip.label), unname(nmAllLbl[1:nTips(p4)]))
    checkEquals(unname(p4c@node.label),
                unname(nmAllLbl[(nTips(p4)+1):(nTips(p4)+nNodes(p4))]))


    p4c <- p4
    tmpNm <- nmAllLbl[order(as.numeric(names(nmAllLbl)))]
    labels(p4c, "all", use.names=TRUE) <- nmAllLbl
    checkTrue(all(names(p4c@tip.label) %in% nodeId(p4c, "tip")))
    checkTrue(all(names(p4c@node.label) %in% nodeId(p4c, "internal")))
    checkEquals(p4c@tip.label, tmpNm[names(tmpNm) %in% nodeId(p4c, "tip")])
    checkEquals(p4c@node.label, tmpNm[names(tmpNm) %in% nodeId(p4c, "internal")])
}

test.labelsAllPhylo4d <- function() {

    allLbl <- paste("n", 1:(nTips(p4d)+nNodes(p4d)), sep="")
    nmAllLbl <- allLbl
    names(nmAllLbl) <- nodeId(p4d, "all")

    allLbl <- sample(allLbl)
    nmAllLbl <- sample(nmAllLbl)

    p4dc <- p4d
    labels(p4dc, "all") <- allLbl
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    checkTrue(all(names(p4dc@node.label) %in% nodeId(p4dc, "internal")))
    checkEquals(unname(p4dc@tip.label), allLbl[1:nTips(p4d)])
    checkEquals(unname(p4dc@node.label),
                allLbl[(nTips(p4d)+1):(nTips(p4d)+nNodes(p4d))])

    p4dc <- p4d
    labels(p4dc, "all") <- nmAllLbl
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    checkTrue(all(names(p4dc@node.label) %in% nodeId(p4dc, "internal")))
    checkEquals(unname(p4dc@tip.label), unname(nmAllLbl[1:nTips(p4d)]))
    checkEquals(unname(p4dc@node.label),
                unname(nmAllLbl[(nTips(p4d)+1):(nTips(p4d)+nNodes(p4d))]))


    p4dc <- p4d
    tmpNm <- nmAllLbl[order(as.numeric(names(nmAllLbl)))]
    labels(p4dc, "all", use.names=TRUE) <- nmAllLbl
    checkTrue(all(names(p4dc@tip.label) %in% nodeId(p4dc, "tip")))
    checkTrue(all(names(p4dc@node.label) %in% nodeId(p4dc, "internal")))
    checkEquals(p4dc@tip.label, tmpNm[names(tmpNm) %in% nodeId(p4dc, "tip")])
    checkEquals(p4dc@node.label, tmpNm[names(tmpNm) %in% nodeId(p4dc, "internal")])
}
