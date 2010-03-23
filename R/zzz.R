# only used when there is no namespace
#.First.lib <- function(lib, pkg) {
#	require(ape)
#       library.dynam("phylobase", pkg, lib )
#}

".phylobase.Options" <-
    list(retic = "warn",
         singleton = "warn",
         multiroot = "warn",
         poly = "warn",
         allow.duplicated.labels = "fail")

# use this with a namespace
.onLoad <- function(lib, pkg) {
    require(ape)
    require(methods)
}

.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".phylobase.Options", asNamespace("phylobase"))
}


