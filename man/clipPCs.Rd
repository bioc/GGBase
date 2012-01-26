\name{clipPCs}
\alias{clipPCs}
\alias{regressOut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
transformations of expression data in smlSet instances
}
\description{
transformations of expression data in smlSet instances
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clipPCs(smlSet, inds2drop, center = TRUE)

regressOut(sms, rhs, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
instance of \code{\linkS4class{smlSet}}
}
  \item{sms}{
instance of \code{\linkS4class{smlSet}}
}
  \item{inds2drop}{
vector of PCs to be eliminated by setting the associated diagonal elements in the
SVD to zero before recomposing the matrix of expression values
}
  \item{center}{
logical, passed to \code{\link[stats]{prcomp}}
}
  \item{rhs}{formula fragment (no dependent variable) used to form residuals
in a reexpression of the expression matrix}
  \item{\dots}{arguments passed to \code{\link[limma]{lmFit}}}
}
\details{
\code{clipPCs} is an operation on the n x p transposed matrix X of expression data.
The singular value decomposition X = UDV^t is formed, the diagonal elements of D
corresponding to \code{inds2drop} are set to zero yielding the diagonal matrix E, and
then Y = UEV^t is computed and transposed to replace the expression data.

\code{regressOut} obtains residuals after genewise regression of expression on
the design matrix specified by the \code{rhs}; \code{\link[limma]{lmFit}} is used to
compute coefficients, linear predictions and residuals.
}
\value{
an instance of \code{\linkS4class{smlSet}}
}

\references{
The use of PCA-based adjustments to remove mass extraneous effects from
expression matrices has been criticized in work of Oliver Stegle and Jeffrey Leek,
who offer Bayesian PEER and SVA respectively as alternative solutions.
The PCA-based method seems to have reasonable 
effectiveness in examples worked with GGdata.
}
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
\dontrun{  # this would induce cyclic dependency, but should
           # run manually
 if ("GGtools" \%in\% installed.packages()[,1]) {
  require("GGtools")
  c20 = getSS("GGtools", "20")
  t1 = gwSnpTests(genesym("CPNE1")~male, c20, chrnum("20"))
  topSnps(t1)
  c20c = clipPCs(c20, 1:10)
  t2 = gwSnpTests(genesym("CPNE1")~male, c20c, chrnum("20"))
  topSnps(t2)
  }
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }