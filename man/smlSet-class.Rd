\name{smlSet-class}
\docType{class}
\alias{smlSet-class}
\alias{smList,smlSet-method}
\alias{smlEnv,smlSet-method}
\alias{smList}
\alias{smlEnv}

\title{Documentation on
  S4 class "smlSet" an eSet-derived container for snpMatrix lists,
 allowing efficient combination of SNP chip genotyping with microarray
 expression data, and allied classes }
\description{Documentation on
  S4 class "smlSet" an eSet-derived container for snpMatrix lists,
 allowing efficient combination of SNP chip genotyping with microarray
 expression data, and allied classes }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("smlSet", assayData, phenoData, featureData, experimentData, annotation, ...)}.
 These objects respond to interrogation on samples, expression values,
   SNP values, and other metadata.
}
\section{Slots}{
	 \describe{
    \item{\code{smlEnv}:}{Object of class \code{"environment"} an environment
with single key \code{smList} pointing to a list of
package snpMatrix snp.matrix instances }
    \item{\code{snpLocPathMaker}:}{Object of class \code{"function"} obsolete }
    \item{\code{chromInds}:}{Object of class \code{"numeric"} numeric
vector indicating what chromosomes are represented in the \code{smlEnv} }
    \item{\code{organism}:}{Object of class \code{"character"} informal,
"Hs" recommended for human }
    \item{\code{snpLocPackage}:}{Object of class \code{"character"} tells
what package houses the netCDF store for SNP location metadata }
    \item{\code{snpLocRef}:}{Object of class \code{"character"} tells
the symbol exported by the snpLocPackage with which the netCDF store connection
can be retrieved.  That symbol refers to an instance of the S3 class
'ncdf', in package 'ncdf'. }
    \item{\code{activeSnpInds}:}{Object of class \code{"numeric"} redundant with
chromInds slot }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} intended
to hold expression data coordinated with the smlEnv data }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} 
   standard sample-level data container from eSet design }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} 
 standard feature-level metadata container, implied usage is for documenting
the expresion data elements }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} standard
metadata container from Biobase eSet design }
    \item{\code{annotation}:}{Object of class \code{"character"} vector
of strings; element exprs gives the Bioconductor annotation package
(.db type) for decoding expression feature identifiers. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} class
version tracking metadata }
  }
}
\section{Extends}{
Class \code{"\link[Biobase:class.eSet]{eSet}"}, directly.
Class \code{"\link[Biobase:class.VersionedBiobase]{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\link[Biobase:class.Versioned]{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{smList}{\code{signature(x = "smlSet")}: retrieves the
actual list of snp.matrix entities }
    \item{smlEnv}{\code{signature(x = "smlSet")}: retrieves the
environment holding snp.matrix entities }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
GGtools package makes extensive use of these classes and methods.
}
\examples{
showClass("smlSet")
}
\keyword{classes}
