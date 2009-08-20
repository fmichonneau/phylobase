# only used when there is no namespace
#.First.lib <- function(lib, pkg) {
#	require(ape)
#       library.dynam("phylobase", pkg, lib )
#}

# use this with a namespace
.onLoad <- function(lib, pkg) {
    require(ape)
    require(methods)
}
