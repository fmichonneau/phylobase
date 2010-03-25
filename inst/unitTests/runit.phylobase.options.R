
###
###  phylobase.options
###

op <- phylobase.options()
test.phylobase.options <- function() {

    ## test match.arg
    checkException(phylobase.options(retic="test"))
    no <- phylobase.options(retic="f")
    checkIdentical(no$retic, "fail")
    phylobase.options(op)

    ## test multiple args
    no <- phylobase.options(retic="f", poly="f")
    checkIdentical(no$retic, "fail")
    checkIdentical(no$poly, "fail")
    phylobase.options(op)

    ## check some failures
    checkException(phylobase.options(1))
    checkException(phylobase.options("foobar"="foo"))
}
