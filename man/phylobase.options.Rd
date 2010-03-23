\name{phylobase.options}
\alias{phylobase.options}
\title{Set or return options of phylobase}
\description{
  Provides a mean to control the validity of \code{phylobase} objects
  such as singletons, reticulated trees, polytomies, etc. 
}
\usage{
phylobase.options(\dots)
}
\arguments{
  \item{\dots}{a list may be given as the only argument, or any number
    of arguments may be in the \code{name=value} form, or no argument at
    all may be given.  See the Value and Details sections for
    explanation.}
  Arguments which are set by a function call will remain in effect until
  the end of the current session, unless overwritten by a
  subsequent call.  In addition, they can be added as optional
  parameters of calls to specific functions of \code{phylobase}; in
  this case, their effect is limited to that function call.  See the
  documentation of specific functions for the list of options which are
  recognised by that function.
  \itemize{
    \item{retic}{\dQuote{warn}, \dQuote{fail} or \dQuote{OK}. Are
      reticulated trees allowed? \dQuote{warn} returns a warning
      (default) and \dQuote{fail} returns an error message}
    
    \item{singleton}{\dQuote{warn}, \dQuote{fail} or \dQuote{OK}. Are
      singleton nodes allowed? \dQuote{warn} returns a warning
      (default) and \dQuote{fail} returns an error message}
    
    \item{multiroot}{\dQuote{warn}, \dQuote{fail} or \dQuote{OK}. Are
      multiple roots allowed? \dQuote{warn} returns a warning
      (default) and \dQuote{fail} returns an error message}
    
    \item{poly}{\dQuote{warn}, \dQuote{fail} or \dQuote{OK}. Are
      polytomies allowed? \dQuote{warn} returns a warning
      (default) and \dQuote{fail} returns an error message}
    
    \item{allow.duplicated.labels}{\dQuote{warn}, \dQuote{fail} or
      \dQuote{OK}. Are duplicated labels allowed? \dQuote{warn} returns
      a warning and \dQuote{fail} (default) returns an error message}
  }
}
\value{
  A list with the updated values of the parameters. If arguments are
      provided, the returned list is invisible.
}
\details{
  The parameter values set via a call to this function will remain in
  effect for the rest of the session, affecting the subsequent behavior
  of phylobase.
}

\examples{
\dontrun{
phylobase.options(poly="fail")
# subsequent trees with polytomies will fail the validity check
}
}
\author{Francois Michonneau (adapted from the package \code{sm})}
\keyword{phylobase}
\keyword{validator}