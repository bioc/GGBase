\name{gwSnpScreenResult-class}
\docType{class}
\alias{gwSnpScreenResult-class}
\alias{cwSnpScreenResult-class}
\alias{multiGwSnpScreenResult-class}
\alias{filteredGwSnpScreenResult-class}
\alias{filteredMultiGwSnpScreenResult-class}
\alias{show,filteredMultiGwSnpScreenResult-method}
\alias{show,filteredGwSnpScreenResult-method}
\alias{show,cwSnpScreenResult-method}
\alias{show,gwSnpScreenResult-method}
\alias{show,multiGwSnpScreenResult-method}
\alias{plot,cwSnpScreenResult,ANY-method}
\alias{plot,gwSnpScreenResult,ANY-method}
\alias{plot,multiGwSnpScreenResult,ANY-method}
\alias{[,cwSnpScreenResult,ANY,ANY,ANY-method}
%\alias{plot,filteredGwSnpScreenResult,ANY-method}
%\alias{plot,filteredMultiGwSnpScreenResult,ANY-method}


\title{Class "gwSnpScreenResult" -- containers for GGtools gwSnpScreen
method outputs and allied objects}
\description{Class "gwSnpScreenResult" -- container for GGtools gwSnpScreen
method outputs and allied objects}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gwSnpScreenResult", ...)}.
These will be primarily lists of inference tables (snps are rows,
columns are statistics and p-values).  Additional slots
manage analysis metadata.

gwSnpScreenResult is intended for genome-wide analysis of expression
for a single gene.

cwSnpScreenResult is intended for the restriction to a single
chromosome.

multiGwSnpScreenResult is intended for analyses with multiple genes.

Because the vast majority of tests are uninformative, early
filtering is important for managing object sizes.
Instances of filteredGwSnpScreenResult and filteredMultiGwSnpScreenResult
are created when a \texttt{snpdepth} parameter is used with
gwSnpTests.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} containing
inference tables (snps are rows, columns are statistics and
p-values) }
    \item{\code{gene}:}{Object of class \code{"character"} typically
the HUGO symbol of the gene analyzed }
    \item{\code{psid}:}{Object of class \code{"character"} the
feature identifier of the associated microarray }
    \item{\code{annotation}:}{Object of class \code{"character"} vector
of relevant annotation package identifier names }
were used}
    \item{\code{formula}:}{Object of class \code{"formula"} the
formula used to fit the model relating expression to
genotype }
 
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
}
\section{Methods}{
plot and show
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("gwSnpScreenResult")
showClass("cwSnpScreenResult")
}
\keyword{classes}
