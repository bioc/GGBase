\name{racExSet-class}
\docType{class}
\alias{racExSet-class}
\alias{initialize,racExSet-method}

\title{Class "racExSet" -- container for genotype data
in 'rare allele count' form and expression data.  Approaching
obsolescence, see smlSet for more efficient representation.}

\description{Class "racExSet" -- container for genotype data
in 'rare allele count' form and expression data.  Approaching
obsolescence, see smlSet for more efficient representation.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("racExSet", phenoData, experimentData, annotation, exprs, racs, rarebase, SNPalleles)}.
}
\section{Slots}{
	 \describe{
    \item{\code{racAssays}:}{Object of class \code{"AssayData"} S x N
matrix holding rare allele counts }
    \item{\code{rarebase}:}{Object of class \code{"character"} nucleotide
code for rare allele of each SNP in the population assayed }
    \item{\code{SNPalleles}:}{Object of class \code{"character"} 
nucleotide code for polymorphism }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} expression
data container }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} 
sample data container }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} expression feature metadata container }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} experiment
level metadata }
    \item{\code{annotation}:}{Object of class \code{"character"} Bioconductor
annotation package name }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} class
version metadata }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "racExSet")}: ... }
	 }
}
%\references{ ~put references to the literature/web site here ~ }
\author{ VJ Carey <stvjc@channing.harvard.edu> }

\seealso{
murine examples in GGtools continue to use this structure	
}
\examples{
showClass("racExSet")
}
\keyword{classes}
