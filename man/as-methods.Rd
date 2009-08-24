\name{as}
\docType{methods}
\alias{as}
\alias{as-method}
\alias{as,phylo,phylo4-method}
\alias{as,phylo,phylo4d-method}
\alias{as,multiPhylo4,multiPhylo-method}
\alias{as,multiPhylo,multiPhylo4-method}
\alias{as,multiPhylo4d,multiPhylo-method}
\alias{as,phylo4,phylo-method}
\alias{as,phylo4d,phylo-method}
\alias{as,phylo4,data.frame-method}
\alias{as,phylo4d,data.frame-method}
\alias{as,phylo4vcov,phylo4-method}
\alias{as,phylo4,phylo4vcov-method}
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
\alias{coerce,phylo4vcov,phylo4-method}
\alias{coerce,phylo4,phylo4vcov-method}

\title{Converting between phylo4/phylo4d and other phylogenetic tree formats}

\section{Usage}{
\code{as(object, class)}
}

\section{Arguments}{
\describe{
  \item{\code{object}}{a tree of class \code{phylo4}, \code{phylo} or
    \code{phylog}, or tree and data object of class \code{phylo4d}.}
  \item{\code{class}}{the name of the class to which \code{tree} should
    be coerced (e.g., \code{"phylo4"} or \code{"data.frame"}).} 
}
}
\description{
  Translation functions to convert between phylobase objects
  (\code{phylo4} or \code{phylo4d}), and objects used by other
  comparative methods packages in R:  \code{ape} objects (\code{phylo},
  \code{multiPhylo}), \code{ade4} objects (\code{phylog}), and to
  \code{data.frame} respresentation.
}

\section{Methods}{
  
  Coerce from one object class to another using
  \code{as(object, "class")}, where the \code{object} is of the old class
  and the returned object is of the new class \code{"class"}. The
  \code{as} function examines the class of \code{object} and the new
  \code{"class"} specified to choose the proper conversion without
  additional information from the user. Conversions exist for
  combinations:

  \describe{
  
    \item{}{\code{phylobase} to \code{phylobase} formats:}
    \itemize{   
      \item{\code{as(object, "phylo4d")}}{where object is of class
	\linkS4class{phylo4} and returns an object of class
	\linkS4class{phylo4d}, with empty data.} 
      \item{\code{as(object, "phylo4")}}{where object is of class
	\linkS4class{phylo4d} and returns an object of class
	\linkS4class{phylo4}. If data are dropped during the conversion, a
	warning message is produced. A similar conversion can be done by
	using the function \code{extractTree}, but in this case, no error
	message is produced.}
    }
    
    \item{}{\code{phylobase} to \code{ape} formats:}
    \itemize{
      \item{\code{as(object, "phylo")}}{where object is of class
	\linkS4class{phylo4} or \linkS4class{phylo4d} and returns an
	object of class \code{phylo}. If data are dropped during the
	conversion from a \code{phylo4d} object, a warning message is
	produced.}
      \item{\code{as(object, "multiPhylo")}}{ ~~Not implemented
	yet. where object is of class \linkS4class{multiPhylo4} and returns
	an object of class \code{multiPhylo}. }
    }
    
    \item{}{\code{ape} to \code{phylobase} formats:}
    \itemize{
      \item{\code{as(object, "phylo4")}}{where object is of class
	\code{phylo} and returns an object of class \code{phylo4}.}
      \item{\code{as(object, "phylo4d")}}{where object is of class
	\code{phylo} and returns an object of class \code{phylo4d}, with
	empty data.} 
      \item{\code{as(object, "multiPhylo4")}}{ ~~Not implemented
	yet. where object is of class \code{multiPhylo} and returns an
	object of class \code{multiPhylo4}. } 
    }
    
    \item{}{\code{phylobase} to \code{ade4} formats:}
    \itemize{
      \item{\code{as(object, "phylog")}}{where object is of class
	\code{phylo4} and returns an object of class \linkS4class{phylog}.} 
    }
    
    \item{}{\code{phylobase} format to \code{data.frame}:}
    \itemize{
      \item{\code{as(object, "data.frame")}}{where object is of class
	\code{phylo4} or \code{phylo4d} and returns an object of class
	\code{data.frame}, with data included in the case of
	\code{phylo4d}.} 
    }
  }
}


\author{Ben Bolker, Thibaut Jombart, Marguerite Butler, Steve Kembel}
\seealso{ generic \code{\link[methods]{as}}, \code{\link{phylo4}},
  \code{\link{phylo4d}}, \code{\link{extractTree}}, the original
  \code{\link[ade4]{phylog}} from the \code{ade4} package and
  \code{\link[ape]{as.phylo}} from the \code{ape} package. 
}


\examples{
tree.owls <- read.tree(text="(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);")
## round trip conversion
tree_in_phylo <- tree.owls                  # tree is a phylo object
(tree_in_phylo4 <- as(tree.owls,"phylo4"))  # phylo converted to phylo4
identical(tree_in_phylo,as(tree_in_phylo4,"phylo"))
## test if phylo, and phlyo4 converted to phylo are identical
## (no, because of dimnames)

## Conversion to phylog (ade4)
as(tree_in_phylo4, "phylog")

## Conversion to data.frame
as(tree_in_phylo4, "data.frame")

## Conversion to phylo (ape) 
as(tree_in_phylo4, "phylo")

## Conversion to phylo4d, (data slots empty)    
as(tree_in_phylo4, "phylo4d")
}
\keyword{methods}
\concept{phylo4 tree formats}
\concept{convert tree formats}
\concept{conversion between tree formats}
