\name{MAFfilter}
\Rdversion{1.1}
\alias{MAFfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
restrict SNP in an smlSet to range of minor allele frequencies (MAF)
}
\description{
restrict SNP in an smlSet to range of minor allele frequencies (MAF)
}
\usage{
MAFfilter(x, lower = 0, upper = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
smlSet instance
}
  \item{lower}{
numeric lower bound on minor allele frequency for keeping a SNP
}
  \item{upper}{
numeric upper bound on minor allele frequency for keeping a SNP
}
}
\details{
uses \code{\link[snpMatrix]{snp.matrix-class}} summary method from snpMatrix
}
\value{
revised instance of \code{\link{smlSet-class}}
}

%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(smlSet.example)
sapply(smList(MAFfilter(smlSet.example, lower=.1)), dim)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
