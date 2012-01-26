\name{multiCisTestResult-class}
\Rdversion{1.1}
\docType{class}
\alias{multiCisTestResult-class}
\alias{show,multiCisTestResult-method}

\title{Class "multiCisTestResult"  }
\description{	object to contain results of restricted gene-centric searches for eQTL } 
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("multiCisTestResult", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} -- list of
results of \code{\link[snpStats]{snp.rhs.tests}} }
    \item{\code{conditions}:}{Object of class \code{"list"} -- list
of runtime conditions encountered }
    \item{\code{call}:}{Object of class \code{"call"} -- for auditing,
the call used is saved }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "multiCisTestResult")}: ... }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("multiCisTestResult")
}
\keyword{classes}
