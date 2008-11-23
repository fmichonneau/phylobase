\name{pdata-class}
\docType{class}
\alias{ptypes}
\alias{pdata-class}
\alias{[<-,pdata-method}
\alias{[,pdata-method}
\alias{[,pdata,ANY,ANY,ANY-method}
\alias{[[,pdata-method}
\alias{[[<-,pdata-method}
\alias{[[,pdata,ANY,ANY-method}
\alias{[[,pdata,ANY,missing-method}
\alias{$,pdata-method}
\alias{$<-,pdata-method}

\title{Class "pdata"}
\description{Data class for phylo4d objects}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("pdata", ...)}.
}
\section{Slots}{
	 \describe{
	   \item{\code{data}:}{A data frame of tip or node data.  Can be
	     accessed transparently with any of the data frame accessor methods}
	   \item{\code{type}:}{A factor with length equal to \code{ncol(data)}
	     and levels ("multitype","binary","continuous","DNA","RNA","aacid")}
	   \item{\code{comment}:}{A character vector of length \code{ncol(data)}}
	   \item{\code{metadata}:}{An arbitrary list, for storing other user-defined metadata}
}
}
\section{Methods}{
  \describe{
    \item{\[}{\code{signature(x = "pdata")}: access data rows, columns or elements}
    \item{\[<-}{\code{signature(x = "pdata")}: set data rows, columns or elements} }
    \item{\$}{\code{signature(x = "pdata")}: access data rows, columns or elements}
    \item{\[\$<-}{\code{signature(x = "pdata")}: set data rows, columns or elements}
    \item{\[\[<-}{\code{signature(x = "pdata")}: set data columns or elements}
    \item{\[\[}{\code{signature(x = "pdata", i = "ANY", j = "ANY")}: access data columns or elements }
    \item{\[\[}{\code{signature(x = "pdata", i = "ANY", j = "missing")}: set data columns or elements }}
  \author{Ben Bolker}
  % commented to eliminate Note: on build
  %\examples{
%}
\keyword{classes}
