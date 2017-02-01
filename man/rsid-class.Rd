\name{rsid-class}
\Rdversion{1.1}
\docType{class}
\alias{rsid-class}
\alias{probeId-class}
\alias{genesym-class}
\alias{rsid,character-method}
\alias{rsid,numeric-method}
\alias{probeId,character-method}
\alias{genesym,character-method}
\alias{genesym}
\alias{rsid}
\alias{probeId}
\alias{casting}

\title{Class \code{"rsid"}}
\description{
rsid(), probeId(), and genesym() are basic casting
methods that assign a type to a token.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rsid", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} that holds the content
to which a type is associated by the method. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
Class \code{"\linkS4class{data.frameRowLabels}"}, by class "character", distance 2.
Class \code{"\linkS4class{SuperClassMethod}"}, by class "character", distance 2.
Class \code{"\linkS4class{characterORconnection}"}, by class "character", distance 2.
Class \code{"\linkS4class{characterORMIAME}"}, by class "character", distance 2.
Class \code{"\linkS4class{atomicVector}"}, by class "character", distance 2.
Class \code{"\linkS4class{index}"}, by class "character", distance 2.
%Class \code{"\linkS4class{atomic}"}, by class "character", distance 2.
%Class \code{"\linkS4class{character_OR_NULL}"}, by class "character", distance 2.
%Class \code{"\linkS4class{vector_OR_factor}"}, by class "character", distance 3.
}
\section{Methods}{
No methods defined with class "probeId" or "genesym" in the signature.
However \code{[} with a \code{\link{smlSet-class}} argument will
dispatch differently if instances of these classes are supplied.
Various GGtools reporting functions will use "rsid" for dispatch.
}
%\references{
%%  ~~put references to the literature/web site here~~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
%}
\examples{
showClass("rsid")
}
\keyword{classes}
