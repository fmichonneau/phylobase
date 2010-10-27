\name{phylog-class}
\docType{class}
\alias{phylog-class}
\alias{coerce,phylo4,phylog-method}
\title{Class "phylog"}
\description{S4 version of the class \code{phylog} from the \code{ade4} package.}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Extends}{
Class \code{\link[ade4]{phylog}}, directly.
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "phylo4", to = "phylog")}:
      a \linkS4class{phylo4} object can be coerced to a \code{phylog} object using \code{as(object,"phylog")}. }
	 }
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\seealso{The original \code{\link[ade4]{phylog}} from the
  \code{ade4} package. 
}
\section{Warning}{\code{phylog} objects from \code{ade4} are now deprecated;
users are encouraged to use the \code{adephylo} package, which uses
\code{phylo} and \code{phylo4} objects directly.}
\keyword{classes}
