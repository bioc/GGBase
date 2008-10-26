\name{smlSet-class}
\docType{class}
\alias{smlSet-class}
\alias{smList,smlSet-method}
\alias{getAlleles,smlSet,rsid-method}
\alias{getAlleles}
\alias{smlEnv,smlSet-method}
\alias{[,snp.matrix,ANY,rsid,ANY-method}
\alias{[,smlSet,ANY,ANY,ANY-method}
\alias{exprs,smlSet-method}
\alias{show,smlSet-method}
\alias{combine,smlSet,smlSet-method}
\alias{snpNames}
\alias{snpNames,smlSet,chrnum-method}
\alias{snpNames,smlSet,missing-method}
%Undocumented code objects:


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
    \item{\code{organism}:}{Object of class \code{"character"} informal,
"Hs" recommended for human }
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
giving the Bioconductor annotation package
(.db type) for decoding expression feature identifiers. }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} class
version tracking metadata }
  }
}
\section{Extends}{
Class \code{\link[Biobase:class.eSet]{eSet}}, directly.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, 
by class "eSet", distance 2.
Class \code{\link[Biobase:class.Versioned]{Versioned}}, by class "eSet", 
distance 3.
}
\section{Methods}{
  \describe{
    \item{smList}{\code{signature(x = "smlSet")}: retrieves the
actual list of snp.matrix entities }
    \item{smlEnv}{\code{signature(x = "smlSet")}: retrieves the
environment holding snp.matrix entities }
    \item{exprs}{\code{signature(x = "smlSet")}: retrieves the
matrix of expression values }
    \item{\code{combine}:}{concatenates expression
data and forms intersection of SNP sets}
    \item{\code{getAlleles(smlSet, rsid)}:}{returns A/B notations for SNP
determined by rsid}
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{ We have included a \code{[} method for snp.matrix instances
that accepts an rsid instance as a column selector. }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
GGtools package makes extensive use of these classes and methods.
}
\examples{
showClass("smlSet")
}
\keyword{classes}
