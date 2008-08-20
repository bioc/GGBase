\name{genesym-class}
\docType{class}
\alias{genesym-class}
\alias{genesym,character-method}
\alias{chrnum-class}
\alias{chrnum,numeric-method}
\alias{chrnum,character-method}
\alias{rsid-class}
\alias{rsid,character-method}
\alias{rsid,numeric-method}
\alias{probeId,character-method}
\alias{cnumOrMissing-class}
\alias{probeId-class}
\alias{genesym}
\alias{probeId}
\alias{chrnum}
\alias{rsid}
\alias{show,rsid-method}
\alias{show,chrnum-method}
\alias{snpdepth}
\alias{snpdepth-class}
\alias{cnumOrMissing}


\title{Class "genesym" and other casting classes}
\description{	 classes that help establish symbol semantics for
dispatching }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("genesym", ...)}.
and generally just extend character or numeric so that vector operations
are straightforward

\code{snpdepth} identifies a number that will be used as the number
of chromosome-specific test results to be retained in any
genome-wide screen
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
Class \code{"\linkS4class{characterORMIAME}"}, by class "character", distance 2.
}
%\section{Methods}{
%No methods defined with class "genesym" in the signature.
%}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\examples{
showClass("genesym")
genesym("CPNE1")
}
\keyword{classes}
