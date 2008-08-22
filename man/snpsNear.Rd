\name{snpsNear}
\alias{snpsNear}
\title{ obtain list of rs numbers for snps near a gene }
\description{obtain list of rs numbers for snps near a gene
}
\usage{
snpsNear(gsym, radius=1e+05)
}
\arguments{
  \item{gsym}{ instance of genesym class [e.g., use genesym(string) for gene 'string']}
  \item{radius}{ number of base-pairs in each direction to look}
}
\details{
simple arithmetic based on output of snpLocs.Hs
}
\value{
character vector of rsxxxxxx, dbSNP id, according to locations from
SNPlocs.Hsapiens.dbSNP.20071016 package
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
\note{first invocation can take longer than subsequent, if snpLocs.Hs has
not been invoked previously  }


%\seealso{  }

\examples{
nearc = snpsNear(genesym("CPNE1"), 10000)
}
\keyword{ models }
