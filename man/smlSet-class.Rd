\name{smlSet-class}
\docType{class}
\alias{smlSet}
\alias{smlSet.example}
\alias{smlSet-class}
\alias{smList,smlSet-method}
\alias{getAlleles,smlSet,rsid-method}
\alias{getAlleles}
\alias{smlEnv,smlSet-method}
\alias{[,SnpMatrix,ANY,rsid,ANY-method}
\alias{[,smlSet,ANY,ANY,ANY-method}
\alias{[,smlSet-method}
\alias{exprs,smlSet-method}
\alias{updateObject,smlSet-method}
\alias{snps}
\alias{snps,smlSet,chrnum-method}
\alias{show,smlSet-method}
\alias{coerce,smlSet,ExpressionSet-method}
\alias{combine,smlSet,smlSet-method}
\alias{snpNames}
\alias{snpNames,smlSet,chrnum-method}
\alias{snpNames,smlSet,missing-method}
%Undocumented code objects:


\alias{smList}
\alias{smlEnv}

\title{Documentation on
  S4 class "smlSet" an eSet-derived container for SnpMatrix lists,
 allowing efficient combination of SNP chip genotyping with microarray
 expression data, and allied classes }
\description{Documentation on
  S4 class "smlSet" an eSet-derived container for SnpMatrix lists,
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
package snpStats SnpMatrix instances }
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
actual list of SnpMatrix entities }
    \item{smlEnv}{\code{signature(x = "smlSet")}: retrieves the
environment holding SnpMatrix entities }
    \item{exprs}{\code{signature(x = "smlSet")}: retrieves the
matrix of expression values }
    \item{snps}{\code{signature(x = "smlSet", chr = "chrnum")}: retrieves the
raw matrix of genotype values (SnpMatrix instance from
snpStats package }
    \item{\code{combine}:}{concatenates expression
data and forms intersection of SNP sets}
    \item{\code{getAlleles(smlSet, rsid)}:}{returns A/B notations for SNP
determined by rsid}
    \item{\code{coerce}:}{extracts exprs, phenoData and annotation and
constructs ExpressionSet}



   \item{[}{\code{signature(x = "smlSet", i = "ANY", j = "ANY", drop = "ANY")}:
Quick methods for subsetting elements of smlSets have been provided.  

If X is an \code{smlSet} instance and G is a vector of class \code{\link{probeId-class}},
then \code{X[G,]} will reduce the expression data to the probes specified
in \code{G}.

If X is an \code{smlSet} instance and G is a vector of class \code{\link{chrnum-class}},
then \code{X[G,]} will reduce the SNP genotype data to the SNPs resident
on chromosomes enumerated
in \code{G}.

If X is an \code{smlSet} instance and G is a vector of class \code{\link{rsid-class}},
then \code{X[G,]} will reduce the SNP genotype data to the SNPs enumerated
in the dbSNP id 
in \code{G}.
 }
}


}
%\references{ ~put references to the literature/web site here ~ }

\author{ VJ Carey <stvjc@channing.harvard.edu> }
\note{ We have included a \code{[} method for SnpMatrix instances
that accepts an rsid instance as a column selector. }

% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
GGtools package makes extensive use of these classes and methods.
}
\examples{
showClass("smlSet")
data(smlSet.example)
smlSet.example
validObject(smlSet.example)
# workout on expression components
dim(exprs(smlSet.example))
fn = featureNames(smlSet.example)[1:10]
fn
ss2 = smlSet.example[ probeId(fn), ] # restrict exprs to set of probes
dim(exprs(ss2))
# workout on SNP components
smList(smlSet.example)
dim(smList(ss2)[[1]])
ss2[ chrnum(21), ]  # trivial restriction of SNP to a chromosome
sn = colnames(smList(ss2)[[1]])[1:20]  # get some dbSNP ids
ss3 = ss2[ rsid(sn), ]   # subset the snps
dim(smList(ss3)[[1]])
dim(smList(ss3)[["21"]]) # check names
ss3
as(snps(ss3, chrnum(21)), "character")[1:5,1:5] # generic codes
as(snps(ss3, chrnum(21)), "numeric")[1:5,1:5]  # number copies of B
as(snps(ss3, chrnum(21)), "matrix")[1:5,1:5]  # raw
}
\keyword{classes}
