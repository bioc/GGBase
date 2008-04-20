\name{snpScreenResult-class}
\docType{class}
\alias{snpScreenResult-class}
\alias{twSnpScreenResult-class}

\title{Class "snpScreenResult" -- container
approaching obsolescence, retained for legacy purposes,
that holds expression vs genotype analysis results; see
gwSnpScreenResult for more current approach. }
\description{Class "snpScreenResult" -- container
approaching obsolescence, retained for legacy purposes,
that holds expression vs genotype analysis results; see
gwSnpScreenResult for more current approach. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("snpScreenResult", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} list of
inference tables }
    \item{\code{call}:}{Object of class \code{"call"} audit information }
    \item{\code{gene}:}{Object of class \code{"character"} gene name }
    \item{\code{locs}:}{Object of class \code{"numeric"} locations
of SNPs }
    \item{\code{chr}:}{Object of class \code{"character"} chromosome
in play }
    \item{\code{fitter}:}{Object of class \code{"ANY"} audit information }
    \item{\code{annotation}:}{Object of class \code{"character"}  propagated
  from expression data }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
No methods defined with class "snpScreenResult" in the signature in
GGBase, but see GGtools.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{ twSnpScreenResult is provided for multiple gene analyses.  See
multiGwSnpScreenResult for the more current approach to such
analyses }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("snpScreenResult")
}
\keyword{classes}
