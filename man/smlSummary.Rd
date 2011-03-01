\name{smlSummary}
\Rdversion{1.1}
\alias{smlSummary}
\alias{smlSummary-class}
\alias{show,smlSummary-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
class and function to summarize frequency information
on genotypes in an smlSet
}
\description{
generates information on sample size, minor allele frequency, specific
call frequencies, and HWE test results on all SNP in an smlSet
}
\usage{
smlSummary(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
instance of \code{\link{smlSet-class}}
}
}
\details{
to control volume of printout a simple list extending class
is defined for show method
}
\value{
Instance of smlSummary class, which simply extends list.
Each list element is a matrix of results provided by
\code{\link[snpStats:SnpMatrix-class]{summary,SnpMatrix-method}}.
}

%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJCarey <stvjc@channing.harvard.edu>
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
smlSummary(smlSet.example)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
