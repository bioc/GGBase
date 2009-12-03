\name{make_smlSet}
\Rdversion{1.1}
\alias{make_smlSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create an smlSet instance from components
}
\description{
create an smlSet instance from components
}
\usage{
make_smlSet(es, sml, organism = "Homo sapiens")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
ExpressionSet instance
}
  \item{sml}{
list of snp.matrix instances
}
  \item{organism}{
string naming organism
}
}
\details{
combines snp.matrix instances with expression data
}
\value{
instance of smlSet class
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(smlSet.example) # here we just show the mechanics from a working smlSet
es = as(smlSet.example, "ExpressionSet")
sl = smList(smlSet.example)
mm = make_smlSet(es, sl)
validObject(mm)
mm
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
