\name{plot_EvG}
\alias{plot_EvG}
\alias{plot_EvG,genesym,rsid,smlSet-method}
\alias{plot_EvG,probeId,rsid,smlSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
display the association between expression values and genotypes in
an smlSet instance
}
\description{
display the association between expression values and genotypes in
an smlSet instance
}
\usage{
plot_EvG(gsym, rsid, sms, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsym}{
instance of class \code{\linkS4class{genesym}} or \code{\linkS4class{probeId}},
casting a string that names a gene (which will be looked up using the
\code{annotation} slot of \code{sms}) or a probe which must be resident
on the array underlying the expression content of \code{sms}
}
  \item{rsid}{
instance of class \code{\linkS4class{rsid}} naming a SNP with
genotype values given among the columns of the \code{smList} components of \code{sms}
}
  \item{sms}{
an instance of \code{\linkS4class{smlSet}}
}
  \item{\dots}{
additional parameter to \code{plot}
}
}
\details{
When the genotype is categorical, will use boxplots; when genotype has
been imputed and includes expected allele counts, will use a scatterplot.
}
\value{
a plot is rendered on the current display
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey <stvjc@channing.harvard.edu>
%%  ~~who you are~~
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
  s20 = getSS("GGtools", "20")
  plot_EvG(genesym("CPNE1"), rsid("rs6060535"), s20)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
