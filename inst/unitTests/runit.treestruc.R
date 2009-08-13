#
# --- Test treestruc functions ---
#
 
test.hasPoly <- function() {
    # construct simple polytomy
    owls <- read.tree(text =
        "((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3);")
    owls$edge <- matrix(c(4,4,4,1,2,3), ncol=2)
    owls$Nnode <- 1
    owls$edge.length <- owls$edge.length[-4]
    tr <- as(owls, "phylo4")
    checkTrue(hasPoly(tr))
}

