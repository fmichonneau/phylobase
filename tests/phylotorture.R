## torture-testing phylo4 objects.
require(phylobase)
set.seed(1001)
p1 <- list()
n <- 10
## don't want to slow down R CMD check by doing this every time:
## n <- 10000
for (i in 1:n) {
##    e2 <- c(sample(1:5,replace=FALSE,size=5),sample(6:10,replace=FALSE,size=5))
##    e1 <- sample(6:10,replace=TRUE
    e <- matrix(sample(1:10,replace=TRUE,size=10),ncol=2)
    p1[[i]] <- try(phylo4(e),silent=TRUE)
}
OKvals <- sapply(p1,class)!="try-error"
table(sapply(p1[!OKvals],as.character))

if (any(OKvals)) {
    p2 <- p1[OKvals]
    length(p2)
    has.poly <- sapply(p2,hasPoly)
    has.sing <- sapply(p2,hasSingle)
    has.retic <- sapply(p2,hasRetic)   
    ##
    if (any(has.sing)) {
        p4 <- p2[has.sing]
        plot(p4[[1]])  ## gives descriptive error
        t2 <- try(plot(collapse.singles(as(p2[[1]],"phylo"))))
        ## "incorrect number of dimensions"
    }
    if (any(!has.sing)) {
        ## first tree without singles -- HANGS!
        ## don't try the plot in an R session you care about ...
        p3 <- p2[!has.sing]
        ## plot(p2[[13]])
    }
}

## elements 8 and 34 are 
## what SHOULD the rules for trees be?

## (a) reduce node numbers to 1 ... N ?
## (b) check: irreducible, non-cyclic, ... ?

## convert to matrix format for checking?

reduce_nodenums <- function(e) {
    matrix(as.numeric(factor(e)),ncol=2)
}

# make an illegal phylo4 object, does it pass checks?
# a disconnected node:

t1 <- read.tree (text="((a,b), (c,(d, e)));")
plot(t1)

broke1 <- t1
broke1$edge[broke1$edge[,2] ==9, 1] <- 9  # disconnect the node, two subtrees, ((a, b), c)  and (d,e)

as(broke1, "phylo4") -> tree   # makes a phylo4  object with no warning
phylo4(broke1$edge)    # constructor makes a phylo4 object with no warning
## error message comes from ape, not phylo? -- AND
##   error is about singles, not disconnected nodes
print(try(plot(tree), silent=TRUE ))  

# root node value != ntips + 1:

broke2 <- t1
broke2$edge[broke2$edge==6] <- 10

## warning, but no error
## plot(broke2)  ## seems to hang R CMD check??
## generates error, but it's about wrong number of tips, not wrong value at root.
print(try(as(broke2, "phylo4"), silent=TRUE))
## error regarding number of tip labels vs edges and nodes
print(try(phylo4(broke2$edge), silent=TRUE))

# switch root node value (6) with next internal node (7):

broke3 <- broke2
broke3$edge[broke3$edge==7] <- 6
broke3$edge[broke3$edge==10] <- 7

## both of the following now fail with
## "root node is not at position (nTips+1)
try(as(broke3,"phylo4") -> tree3)  # works with no error message
try(phylo4(broke3$edge))    # works with no error message
## plot(tree3)  # would work if we could create it?


# tips have larger numbers than root node:

broke4 <- t1
broke4$edge[broke4$edge==1] <- 11
broke4$edge[broke4$edge==2] <- 12
broke4$edge[broke4$edge==3] <- 13
broke4$edge[broke4$edge==4] <- 14
broke4$edge[broke4$edge==5] <- 15

print(try(as(broke4, "phylo4"), silent=TRUE) )  # error message saying tree has more than one root
print(try(phylo4(broke4$edge),silent=TRUE))     # error message saying tree has more than one root
# print(try(plot(broke4), TRUE))   ## CAUSES R TO HANG!

###
foo <- new('phylo4')
set.seed(1001)
foo@edge <- rcoal(10)$edge
print(try(plot(foo)))

foo@tip.label <- rep('blah',10)
foo@node.label <- rep("",9)

#####
## tree with only 2 tips: will fail under previous versions
## with "Error in if (which(nAncest == 0) != nTips + 1) { : 
##  argument is of length zero"

edge <- matrix(c(3,1,3,2),byrow=TRUE,ncol=2)
p2 <- phylo4(edge)
