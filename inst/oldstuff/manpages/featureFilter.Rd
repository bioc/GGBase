\name{featureFilter}
\Rdversion{1.1}
\alias{featureFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
remove unannotated or undesired features from an smlSet instance
}
\description{
remove unannotated or undesired features from an smlSet instance
}
\usage{
featureFilter(x, requires = c("loc", "autosomal"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
instance of smlSet class
}
  \item{requires}{
character vector -- if "loc" is present, require that a
non-NA value is present in CHRLOC for each feature;
if "autosomal" is present, require that CHR value is in 1:22
(presently assumes human genome)
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
revised smlSet instance excluding features no
}
%references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
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
dim(exprs(smlSet.example))
fff = featureFilter(smlSet.example)
dim(exprs(fff))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
