\name{smlSet-class}
\Rdversion{1.1}
\docType{class}
\alias{smlSet-class}
\alias{[,smlSet,ANY,ANY,ANY-method}
\alias{combine,smlSet,smlSet-method}
\alias{exprs,smlSet-method}
\alias{nsFilter,smlSet-method}
\alias{smList,smlSet-method}
\alias{smList}
\alias{permEx,smlSet-method}
\alias{permEx}

\title{Class \code{"smlSet"}}
\description{
Integrative container for expression plus genotype data.  Genotypes
are stored in an efficient format defined in the snpStats package.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("smlSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
The \code{make\_smlSet} function can also be used
to build smlSet instances.
}
\section{Slots}{
  \describe{
    \item{\code{smlEnv}:}{Object of class \code{"environment"} that has
a key \code{smList} element, to which a list of
\code{\link{SnpMatrix-class}} instances is bound.}
    \item{\code{annotation}:}{Object of class \code{"character"}, describes
featureNames component in terms of the name of the
annotation package that can be used to decode expression probe names. }
    \item{\code{organism}:}{Object of class \code{"character"}, a conventional
string. }
    \item{\code{assayData}:}{Object of class \code{"AssayData"},
manages the expression data.  See \code{\link{AssayData-class}}.  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ,
manages sample level data. See \code{\link{AnnotatedDataFrame-class}}.}
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"},
manages metadata on expression probes.  }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"}, manages
metadata on experiment as a whole.  See \code{\link{MIAxE-class}}.  }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"},
additional storage for experimental protocol description.
See \code{\link{eSet-class}}.  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"};
internal management of class version.  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "smlSet", i = "ANY", j = "ANY", drop = "ANY")}: 
will restrict the content of the smlSet instance according to features
of the arguments supplied.  If \code{x} is numeric or a \code{\link{probeId-class}},
the expression content will be restricted.  If \code{y} is numeric or
character, samples will be restricted accordingly. }
    \item{combine}{\code{signature(x = "smlSet", y = "smlSet")}: This method
attempts to amalgamate two smlSet instances in the appropriate way --
assuming that samples are disjoint. }
    \item{nsFilter}{\code{signature(eset = "smlSet")}: executes genefilter's
method (see \code{\link{nsFilter,ExpressionSet-method}}) on the
expression content, and then propagates the additional genotype and
sample level content unchanged. }
    \item{smList}{\code{signature(x = "smlSet")}: retrieves the list
of SnpMatrix instances defining the genotype content. }
    \item{exprs}{\code{signature(x = "smlSet")}: retrieves the matrix
of expression values }
    \item{permEx}{\code{signature(sms = "smlSet")}: uses \code{sample()}
to generate a permutation of sample indices so that expression data
are permuted against genotype data, but original unpermuted sample
identifiers are preserved; this is necessitated by the rematching
behavior of \code{\link{snp.rhs.tests}}}.
	 }
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
showClass("smlSet")
}
\keyword{classes}
