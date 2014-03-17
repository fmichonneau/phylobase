

#' Utilities and Tools for Phylogenetics
#' 
#' Base package for phylogenetic structures and comparative data.
#' 
#' \itemize{ \itemPackage:phylobase \itemType:Package \itemDate:2009
#' \itemDepends:methods, grid, ape(>= 2.1) \itemSuggests:adephylo, MASS
#' \itemLicense:GPL Version 2 or later \itemAuthors:R Hackathon et al.
#' (alphabetically: Ben Bolker, Marguerite Butler, Peter Cowan, Damien de
#' Vienne, Thibaut Jombart, Steve Kembel, Francois Michonneau, David Orme,
#' Brian O'Meara, Emmanuel Paradis, Jim Regetz, Derrick Zwickl )
#' \itemURL:\url{http://phylobase.r-forge.r-project.org/} }
#' 
#' @name phylobase-package
#' @aliases phylobase-package phylobase
#' @docType package
#' @section MoreInfo: See the help index \code{help(package="phylobase")} and
#' run \code{vignette("phylobase", "phylobase")} for further details and
#' examples about how to use \code{phylobase}.
#' @keywords package
NULL


#' Data from Darwin's finches
#' 
#' Phylogenetic tree and morphological data for Darwin's finches, in different
#' formats
#' 
#' 
#' @name geospiza
#' @aliases geospiza geospiza_raw
#' @docType data
#' @format \code{geospiza} is a \code{phylo4d} object; \code{geospiza_raw} is a
#' list containing \code{tree}, a \code{phylo} object (the tree), \code{data},
#' and a data frame with the data (for showing examples of how to merge tree
#' and data)
#' @note Stolen from Luke Harmon's Geiger package, to avoid unnecessary
#' dependencies
#' @source Dolph Schluter via Luke Harmon
#' @keywords datasets
#' @examples
#' 
#' data(geospiza)
#' plot(geospiza)
#' 
NULL



#' 'Owls' data from ape
#' 
#' A tiny tree, for testing/example purposes, using one of the examples from
#' the \code{ape} package
#' 
#' 
#' @name owls4
#' @docType data
#' @format This is the standard 'owls' tree from the \code{ape} package, in
#' \code{phylo4} format.
#' @source From various examples in the \code{ape} package
#' @keywords datasets
#' @examples
#' 
#' data(owls4)
#' plot(owls4)
#' 
NULL
