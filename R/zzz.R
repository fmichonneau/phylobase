# only used when there is no namespace
#.First.lib <- function(lib, pkg) {
#	require(ape)
#       library.dynam("phylobase", pkg, lib )
#}

".phylobase.Options" <-
    list(retic = "warn",
         singleton = "warn",
         multiroot = "warn",
         poly = "ok",
         allow.duplicated.labels = "warn")

# use this with a namespace
## 12 Nov 2012 obsolete?
## .onLoad <- function(lib, pkg) {
##     ## require(ape)
##     require(methods)
## }

.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".phylobase.Options", asNamespace("phylobase"))
}


