\name{oGtypeExSet-class}
\docType{class}
\alias{oGtypeExSet-class}

\title{Class "oGtypeExSet" oligo-based genotyping output + expression
set container -- obsolete}
\description{Class "oGtypeExSet" oligo-based genotyping output + expression
set container -- obsolete}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("oGtypeExSet", assayData, phenoData, featureData, experimentData, annotation, exprs, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{snpCalls}:}{Object of class \code{"list"} ~~ }
    \item{\code{dbConns}:}{Object of class \code{"list"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{ExpressionSet}"}, directly, with explicit coerce.
Class \code{"\linkS4class{eSet}"}, by class "ExpressionSet", distance 2, with explicit coerce.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "ExpressionSet", distance 3, with explicit coerce.
Class \code{"\linkS4class{Versioned}"}, by class "ExpressionSet", distance 4, with explicit coerce.
}
\section{Methods}{
No methods defined with class "oGtypeExSet" in the signature.
}
\section{Note}{
This class was designed to combine oligo SNP chip outputs with
expression data.  This specific workflow is not supported at this
time.  Convert the oligo outputs to snp.matrix and use smlSet
container instead.
}
\examples{
showClass("oGtypeExSet")
}
\keyword{classes}
