\name{snpLocs.Hsapiens}
\alias{snpLocs.Hsapiens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
acquire SNP locations
}
\description{
acquire SNP locations
}
\usage{
snpLocs.Hsapiens(rsid, chrtok, spack = "SNPlocs.Hsapiens.dbSNP.20090506")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rsid}{
dbSNP identifiers
}
  \item{chrtok}{
token to identify chromosomes for SNPlocs package in use
}
  \item{spack}{
name of a SNPlocs package
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
named vector of locations
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
snpLocs.Hsapiens("rs6060535", "chr20")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
