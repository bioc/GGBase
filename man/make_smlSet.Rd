\name{make_smlSet}
\alias{make_smlSet}
\alias{getSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
construct an smlSet instance from existing resources
}
\description{
construct an smlSet instance from existing resources, either
using ExpressionSet and SnpMatrix instances, or a suitably
structured package
}
\usage{
make_smlSet(es, sml, organism = "Homo sapiens", harmonizeSamples = FALSE)

getSS(packname, chrs, renameChrs, probesToKeep=NULL, exFilter=function(x)x,
   wrapperEndo=NULL, checkValid=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
instance of \code{\link{ExpressionSet-class}}
}
  \item{sml}{
named list of \code{\link{SnpMatrix-class}} instances
}
  \item{organism}{
conventional token for species
}
  \item{harmonizeSamples}{
logical indicating whether steps should be taken to be sure that
the components represent identical sets of samples.  can be time
consuming so defaults to FALSE, in which case the user must be sure
that the genotype and expression components are compatible in terms
of sample content
}
  \item{packname}{string naming the installed package from
which expression and genotype data will be acquired}
  \item{chrs}{character vector naming the prefixes of genotype files
to be used in the resulting smlSet instance; see notes below}
  \item{renameChrs}{character vector of same length as \code{chrs} specifying
one-to-one renaming operation for genotype components in smlSet}
  \item{probesToKeep}{character vector identifying probes to be retained
in the constructed smlSet -- to be deprecated in favor of \code{exFilter}}
  \item{exFilter}{function that should accept and return \code{ExpressionSet},
will be executed just after \code{probesToKeep} filter if present}
  \item{wrapperEndo}{function that accepts and returns an smlSet instance,
allowing any sort of transformation of contents acquired with \code{getSS},
executed after all assembly and filtering completed.}
 \item{checkValid}{logical, if TRUE, function will fail if
   created smlSet instance does not pass validObject()}
}
\details{
Packages that work with \code{getSS} can be created out of existing
\code{smlSet} instances using \code{\link{externalize}}.
}
\value{
Instance of \code{\link{smlSet-class}}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
if ("GGtools" \%in\% installed.packages()[,1]) {
 s20 = getSS("GGtools", "20", renameChrs="chr20")
 s20
 make_smlSet( as(s20, "ExpressionSet"), smList(s20) )
}
}
\keyword{ models }
