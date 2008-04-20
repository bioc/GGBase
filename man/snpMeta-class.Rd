\name{snpMeta-class}
\docType{class}
\alias{snpMeta-class}
\alias{snpMetaWhole-class}

\title{Class "snpMeta" representing SNP locations etc., obsolete
for humans, still used for murine data.}
\description{Class "snpMeta" representing SNP locations etc., obsolete
for humans, still used for murine data.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("snpMeta", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{meta}:}{Object of class \code{"environment"} a data
frame of SNP-oriented metadata is stored in this environment,
with key name meta :( }
    \item{\code{chromosome}:}{Object of class \code{"character"} indicator
of relevant chromosome .  snpMetaWhole, for
murine applications, manages container for whole genome.}
  }
}
\section{Methods}{
No methods defined with class "snpMeta" in the signature in
GGBase.  But see GGtools.
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu>}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("snpMeta")
}
\keyword{classes}
