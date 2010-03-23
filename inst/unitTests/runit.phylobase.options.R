
op <- phylobase.options()
test.phylobase.options <- function() {

    ## test match.arg
    checkException(phylobase.options(retic="test"))
    no <- phylobase.options(retic="f")
    checkIdentical(no$retic, "fail")

    ## test multiple args
    phylobase.options(op)
    no <- phylobase.options(retic="f", poly="f")
    checkIdentical(no$retic, "fail")
    checkIdentical(no$poly, "fail")

    ## check some failures
    checkException(phylobase.options(1))
    checkException(phylobase.options("foobar"="foo"))

}

phylobase.options(op)
