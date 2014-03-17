\name{multiPhylo-class}
\docType{class}
\alias{multiPhylo-class}
\alias{multiPhylo4-class}
\alias{multiPhylo4d-class}
\alias{tbind}
\title{multiPhylo4 and extended classes}
\description{Classes for lists of phylogenetic trees.  These classes and methods are planned for a future version of \code{phylobase}.}

% \section{Objects from the Class}{
%   Objects can be created by calls of the form \code{new("phylo4", ...)}
%   Translation from other phylogenetic packages ...
%   }
% \section{Slots}{
%   \describe{
%     \item{\code{phylo}:}{Matrix of edges}
%     \item{\code{edge.length}:}{Numeric vector of edge lengths}
%     \item{\code{Nnode}:}{Integer: number of nodes}
%     \item{\code{tip.label}:}{Character vector of tip labels}
%     \item{\code{root.edge}:}{Integer: root edge (NA if none)}
%     \item{\code{phylolist}:}{(multiPhylo4) list of phylogenetic trees (each of which
%     may be of type \code{phylo4d}, containing node and/or edge data)}
% }
% }
% \section{Methods}{
%   \describe{
%     \item{\$}{\code{signature(x = "phylo4")}: extract a slot}
%     \item{print}{\code{signature(x = "phylo4")}: print method}
%     \item{show}{\code{signature(object = "phylo4")}: show method }
%     \item{summary}{\code{signature(object = "phylo4")}: summary method}
%     \item{names}{\code{signature(x = "phylo4")}: gives the slots names}
% }
% }
% % commented to eliminate Note: on build
% %\usage{
% %  }
%   \arguments{
%     \item{x}{a phylo4 object}
%     \item{printlen}{Maximum number of tip labels to print}
%     \item{quiet}{a logical stating whether the results of the summary
%       should be printed to the screen (FALSE, default) or not (TRUE)}
%     \item{\dots}{additional arguments (for generic consistency)}
%   }
%   \value{
%     The \code{summary} method invisibly returns a list with the
%     following components: \cr
%     - \code{name}: the name of the object \cr
%     - \code{nb.tips}: the number of tips\cr
%     - \code{nb.nodes}: the number of nodes \cr
%     - \code{mean.el}: mean of edge lengths \cr
%     - \code{var.el}: variance of edge lengths (estimate for population) \cr
%     - \code{sumry.el}: summary (i.e. range and quartiles) of the edge
%     lengths \cr
% 
%    The \code{names} method returns a vector of characters corresponding
%    to the names of the slots.
%   }
%   \author{Ben Bolker, Thibaut Jombart}   
% \examples{
% ### An extract from Sibley and Ahlquist (1990)
% cat("(((Strix_aluco:4.2,Asio_otus:4.2):3.1,Athene_noctua:7.3):6.3,Tyto_alba:13.5);", file = "ex.tre", sep = "\n")
% tree.owls <- ape::read.tree("ex.tre", keep.multi = TRUE)
% tree.owls
% (as(tree.owls, "multiPhylo4"))
% }
\keyword{classes}
