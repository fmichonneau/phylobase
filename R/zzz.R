.First.lib <- function(lib, pkg) {
	require(ape)
       library.dynam("phylobase", pkg, lib )
}


