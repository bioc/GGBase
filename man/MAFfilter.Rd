\name{MAFfilter}
\alias{MAFfilter}
\alias{GTFfilter}
\alias{dropMonomorphies}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter genotype contents of an smlSet according to
certain SNP allele frequency features.
}
\description{
Filter genotype contents of an smlSet according to
certain SNP allele frequency features.
}
\usage{
MAFfilter(x, lower = 0, upper = 1)

GTFfilter(x, lower = 0)

dropMonomorphies(sms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{\link{smlSet-class}} instance
}
  \item{sms}{
\code{\link{smlSet-class}} instance
}
  \item{lower}{
lower bound on MAF or GTF to allow retention of associated locus
}
  \item{upper}{
upper bound on MAF or GTF to allow retention of associated locus
}
}
\details{
uses \code{\link[snpStats]{col.summary}} to compute MAF or GTF.  
\code{dropMonomorphies} also uses col.summary.
}
\value{
\code{\link{smlSet-class}} instance
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
if (file.exists(system.file("parts/20.rda", package="GGtools"))) {
 c20 = getSS("GGtools", "20")
 c20
 c20f = MAFfilter(c20, lower=.05)
 c20f
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
