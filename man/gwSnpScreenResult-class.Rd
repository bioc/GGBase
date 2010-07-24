\name{gwSnpScreenResult-class}
\docType{class}
\alias{gwSnpScreenResult}
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
\alias{plot,cwSnpScreenResult,character-method}
\alias{plot,cwSnpScreenResult,missing-method}
\alias{plot,gwSnpScreenResult,ANY-method}
\alias{plot,multiGwSnpScreenResult,ANY-method}
\alias{[,cwSnpScreenResult,ANY,ANY,ANY-method}
%\alias{plot,filteredGwSnpScreenResult,ANY-method}
%\alias{plot,filteredMultiGwSnpScreenResult,ANY-method}
\alias{[,gwSnpScreenResult,ANY,ANY,ANY-method}
\alias{combine,filteredMultiGwSnpScreenResult,filteredMultiGwSnpScreenResult-method}
\alias{combine,multiGwSnpScreenResult,multiGwSnpScreenResult-method}

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
are created when a snpdepth parameter is used with
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
    \item{\code{formula}:}{Object of class \code{"formula"} the
    formula used to fit the model relating expression to
    genotype }
}
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{\link[Biobase:class.AssayData]{AssayData}}, by class "list", 
 distance 2.
}
\section{Methods}{
plot and show
}
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\examples{
showClass("gwSnpScreenResult")
showClass("cwSnpScreenResult")
}
\keyword{classes}
