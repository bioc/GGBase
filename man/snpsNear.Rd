\name{snpsNear}
\alias{snpsNear}
\title{ obtain list of rs numbers for snps near a gene }
\description{obtain list of rs numbers for snps near a gene
}
\usage{
snpsNear(sym, radius=1e+05, chrnum, ...)
}
\arguments{
  \item{sym}{ instance of genesym class [e.g., use genesym(string) for
    gene 'string'], 
or of rsid class, or of numeric class.  An instance of
\code{\link[GSEABase:GeneSet-class]{GeneSet-class}} can also be supplied if it
has \code{geneIdType} AnnotationIdentifier.}
  \item{radius}{ number of base-pairs in each direction to look}
  \item{chrnum}{ chrnum instance .. optional}
  \item{\dots}{ options not now in use }
 }

\details{
simple arithmetic based on output of snpLocs.Hs
}
\value{
character vector of rsxxxxxx, dbSNP id, according to locations from
SNPlocs.Hsapiens.dbSNP.20071016 package, as transferred to snpLocs.Hs
resource in GGBase

note that an attribute 'target' is returned, a named vector
with components chr and loc describing chromosome and location
of the target for which nearby SNPs are sought
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{first invocation can take longer than subsequent, if snpLocs.Hs has
not been invoked previously  }


%\seealso{  }

\examples{
nearc = snpsNear(genesym("BACH1"), 10000, chrnum(21))
data(smlSet.example)
ss = smList(smlSet.example)[[1]]
# following calculation requires new "[" for j an instance of rsid
clo = ss[ , rsid(snpsNear(rsid("rs6060535"), rad=1500, chrnum(20))) ]
clo
# try a gene set
library(GSEABase)
s1 = GeneSet(c("BACH1", "ATP5O"), geneIdType=SymbolIdentifier())
s2 = s1
geneIdType(s2) = AnnotationIdentifier("illuminaHumanv1.db")
s2
sapply(snpsNear(s2), length)
}
\keyword{ models }
