## torture-testing phylo4 objects.
library(phylobase)
set.seed(1001)
p1 <- list()
n <- 10
## don't want to slow down R CMD check by doing this every time:
## n <- 10000
for (i in 1:n) {
    e <- matrix(sample(1:10,replace=TRUE,size=10),ncol=2)
    p1[[i]] <- try(phylo4(e),silent=TRUE)
}
OKvals <- sapply(p1,class)!="try-error"
table(sapply(p1[!OKvals],as.character))

if (any(OKvals)) {
    p2 <- p1[OKvals]
    length(p2)
    has.poly <- sapply(p2,hasPoly)
    has.sing <- sapply(p2,hasSingles)
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
## what SHOULD the rules for trees be?

## (a) reduce node numbers to 1 ... N ?
## (b) check: irreducible, non-cyclic, ... ?

## convert to matrix format for checking?

reduce_nodenums <- function(e) {
    matrix(as.numeric(factor(e)),ncol=2)
}



