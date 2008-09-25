\name{snpsNear}
\alias{snpsNear}
\title{ obtain list of rs numbers for snps near a gene }
\description{obtain list of rs numbers for snps near a gene
}
\usage{
snpsNear(sym, radius=1e+05, chrnum)
}
\arguments{
  \item{sym}{ instance of genesym class [e.g., use genesym(string) for gene 'string'],
or of rsid class, or of numeric class}
  \item{radius}{ number of base-pairs in each direction to look}
  \item{chrnum}{ chrnum instance .. optional}
}
\details{
simple arithmetic based on output of snpLocs.Hs
}
\value{
character vector of rsxxxxxx, dbSNP id, according to locations from
SNPlocs.Hsapiens.dbSNP.20071016 package, as transferred to snpLocs.Hs
resource in GGBase
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{first invocation can take longer than subsequent, if snpLocs.Hs has
not been invoked previously  }


%\seealso{  }

\examples{
nearc = snpsNear(genesym("CPNE1"), 10000, chrnum(20))
library(GGtools)
data(hmceuB36.2021)
ss = smList(hmceuB36.2021)[[1]]
# following calculation requires new "[" for j an instance of rsid
clo = ss[ , rsid(snpsNear(rsid("rs6060535"), rad=1500, chrnum(20))) ]
clo
}
\keyword{ models }
