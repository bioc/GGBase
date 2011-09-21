\name{getSS}
\alias{getSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
construct a small-footprint smlSet instance from a specially structured package
}
\description{
construct a small-footprint smlSet instance from a specially structured package
}
\usage{
getSS(packname, chrs, renameChrs = NULL, probesToKeep = NULL, wrapperEndo = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{packname}{
string naming a package with eset.rda in data folder, defining ExpressionSet instance ex,
and folder parts (from inst in pre-installed image) containing a collection of rda files
holding snpStats \code{\link[snpStats]{SnpMatrix-class}} instances
}
  \item{chrs}{
vector of strings of names of SnpMatrix instances to be included, typically these
are the basenames of files in the \code{parts} folder, and correspond to chromosomes
}
  \item{renameChrs}{
vector of strings of same length of \code{chrs} that will supply names to the
elements of the \code{smList} component of the returned smlSet
}
  \item{probesToKeep}{
vector of strings of probe names to be retained in the returned smlSet
}
  \item{wrapperEndo}{
function receiving and returning an smlSet instance to be invoked prior to
returning the smlSet
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
an instance of \code{\link{smlSet-class}}
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
hm20 = getSS("GGdata", "20", renameChrs="chr20")
}

### The function is currently defined as
#function (packname, chrs, renameChrs = NULL, probesToKeep = NULL, 
#    wrapperEndo = NULL) 
#{
#    if (!is.null(renameChrs) && (length(chrs) != length(renameChrs))) 
#        stop("renameChrs must have same length as chrs in call to getSS")
#    require(packname, character.only = TRUE)
#    ex = get(load(system.file(package = packname, "data/eset.rda")))
#    if (!is.null(probesToKeep)) 
#        ex = ex[probesToKeep, ]
#    partsfol = system.file("parts", package = packname)
#    chk = sapply(chrs, function(x) file.exists(paste(partsfol, 
#        "/", x, ".rda", sep = "")))
#    if (!all(chk)) {
#        cat("requesting ", paste(chrs, ".rda", sep = ""))
#        cat(" but finding\n")
#        print(dir(partsfol))
#        stop("cannot retrieve requested SNP file.")
#    }
#    sml = lapply(chrs, function(x) get(load(paste(partsfol, "/", 
#        x, ".rda", sep = ""))))
#    if (is.null(renameChrs)) 
#        names(sml) = chrs
#    else names(sml) = renameChrs
#    ans = make_smlSet(ex, sml, harmonizeSamples = TRUE)
#    if (is.null(wrapperEndo)) 
#        return(ans)
#    else {
#        ans = wrapperEndo(ans)
#        if (isTRUE(tst <- validObject(ans))) 
#            return(ans)
#        stop(tst)
#    }
#  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
