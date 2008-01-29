\name{coerce-methods}
\docType{methods}
\alias{coerce-methods}

\alias{coerce,phylo,phylo4-method}
\alias{coerce,phylo,phylo4d-method}
\alias{coerce,multiPhylo4,multiPhylo-method}
\alias{coerce,multiPhylo,multiPhylo4-method}
\alias{coerce,multiPhylo4d,multiPhylo-method}
\alias{coerce,phylo4,phylo-method}
\alias{coerce,phylo4d,phylo-method}
\alias{coerce,phylo4,data.frame-method}
\alias{coerce,phylo4d,data.frame-method}

\title{Converting between phylo4 and other phylogenetic trees}
\description{
  Translation functions to convert between phylobase trees or trees+data objects (\code{phylo4}, \code{phylo4d}), and objects used by other comparative methods packages in R:  \code{ape} objects
  (\code{phylo}, \code{multiPhylo}) \code{ade4} objects (\code{phylog}), and to \code{data.frame} respresentation. 
}
\section{Methods}{
  \item{coerce}{from one object class to another using \code{as(object,"Class")}, where the \code{object} is of the old class and the returned object is of the new class \code{"Class"}. The \code{as} function examines the class of \code{object} and the new \code{"Class"} specified to choose the proper conversion without additional information from the user. Conversions exist for combinations:
  \describe{

    \item{\code{phylobase} to \code{phylobase} formats:}{
      \describe{
        \item{\code{as(object,"phylo4d")}}{where object is of class \linkS4class{phylo4} and returns an object of class \linkS4class{phylo4d}, with empty data.}
        \item{\code{as(object,"phylo4")}}{where object is of class \linkS4class{phylo4d} and returns an object of class \linkS4class{phylo4}, dropping data with warning.}
      }
    }
    \item{\code{phylobase} to \code{ape} formats:}{
      \describe{
        \item{\code{as(object,"phylo")}}{where object is of class \linkS4class{phylo4} or \linkS4class{phylo4d} and returns an object of class \code{phylo}, with data dropped in the case of \code{phylo4d} conversion (with warning).}
        \item{\code{as(object,"multiPhylo")}}{ ~~Not implemented yet. where object is of class \linkS4class{multiPhylo4} and returns an object of class \code{multiPhylo}. }
      }
    }
    \item{\code{ape} to \code{phylobase} formats:}{
      \describe{
        \item{\code{as(object,"phylo4")}}{where object is of class \code{phylo} and returns an object of class \code{phylo4}.}
        \item{\code{as(object,"phylo4d")}}{where object is of class \code{phylo} and returns an object of class \code{phylo4d}, with empty data.}
        \item{\code{as(object,"multiPhylo4")}}{ ~~Not implemented yet. where object is of class \code{multiPhylo} and returns an object of class \code{multiPhylo4}. }
      }
    }
    \item{\code{phylobase} to \code{ade4} formats:}{
      \describe{
        \item{\code{as(object,"phylog")}}{where object is of class \code{phylo4} and returns an object of class \linkS4class{phylog}.}
      }
    }
    \item{\code{phylobase} format to \code{data.frame}:}{
      \describe{
        \item{\code{as(object,"data.frame")}}{where object is of class \code{phylo4} or \code{phylo4d} and returns an object of class \code{data.frame}, with data included in the case of \code{phylo4d}.}
      }
    }
  }}
}

\author{Ben Bolker, Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}, Marguerite Butler}
\seealso{ \code{phylo4}, \code{phylo4d}, \code{extract.tree}, the original \code{\link[pkg:ade4]{phylog}} from the
  \code{ade4} package and \code{\link[pkg:ape]{phylo}} from the
  \code{ape} package. }
\examples{
    library(ape)
  example(read.tree)
## round trip conversion 
  tree_in_phylo = tree.owls               # tree is a phylo object 
  tree_in_phylo4 = as(tree.owls,"phylo4")  # phylo converted to phylo4 
  tree_in_phylo4
  identical(tree_in_phylo,as(tree_in_phylo4,"phylo")) # test if phylo, and phlyo4 converted to phylo are identical

  as(tree_in_phylo4, "phylog")     # conversion to phylog (ade4)
  as(tree_in_phylo4, "data.frame")  # conversion to data.frame
  as(tree_in_phylo4, "phylo")      # conversion to phylo (ape)
  as(tree_in_phylo4, "phylo4d")    # conversion to phylo4d, but without data so data slots empty.
}
\keyword{methods}

