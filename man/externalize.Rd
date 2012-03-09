\name{externalize}
\alias{externalize}
\alias{allsnps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a package with a decomposed smlSet instance from
a unified smlSet instance, to reduce memory footprints
}
\description{
create a package with a decomposed smlSet instance from
a unified smlSet instance, to reduce memory footprints
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
externalize(smlSet, packname, 
   author = "Replace Me <auth@a.b.com>", 
   maintainer = "Replace Me <repl@a.b.com>")

allsnps(packname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{smlSet}{
an smlSet instance to be regarded as a source or template for
a package that can be used with greater efficiency of access to
genotype contents
}
  \item{packname}{
name of the package to be generated
}
  \item{author}{
string to be used in DESCRIPTION file of generated package
}
  \item{maintainer}{
string to be used in DESCRIPTION file of generated package
}
}
\details{
The genotype content of the input smlSet is separated into separate RDA files
in the inst/parts folder of the package to be generated.  The ExpressionSet element
of the input smlSet is stored as object \code{ex}
in file \code{eset.rda} in the data folder of the package
to be generated.
}
\value{
creates folder structure and metadata for an installable R package
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
}
\note{
allsnps() gives the list of vectors of snpnames in use in the package
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
 \dontrun{
 setwd(tempdir())
 sms = getSS("GGtools", c("20", "21"))
 externalize(sms, "demopack")
 dir()
 dir("demopack")
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
