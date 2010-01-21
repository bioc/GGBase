\name{snpLocs.Hs}
\alias{snpLocs.Hs}
\alias{hsSnpLocs}
\alias{snpLocs.Hs}
\alias{getSnpLocs}
\alias{getSnpLocs,smlSet,chrnum-method}
\alias{getSnpLocs,smlSet,missing-method}
\alias{snpLocs.Hs,chrnum.rsid-method}
\alias{snpLocs.Hs,chrnum.missing-method}
\alias{snpLocs.Hs,missing,rsid-method}
\alias{snpLocs.Hs,rsid-method}
\alias{snpLocs.Hs,rsid,ANY-method}
\alias{snpLocs-class}
\alias{show,snpLocs-method}
\alias{snpLocs.Hs,chrnum,rsid-method}
\alias{snpLocs.Hs,chrnum,missing-method}
\alias{snpLocs.Hs,missing,rsid-method}
\alias{snpLocs.Hs,sid,ANY-method}

\title{ SNP location accessor }
\description{SNP location accessor
}
\usage{
snpLocs.Hs(cnum, rsid)
}
\arguments{
  \item{cnum}{ chrnum instance indicating chromosome}
  \item{rsid}{ rsid instance giving dbSNP ids for snps of interest}
}
\details{
The SNPlocs.Hsapiens.dbSNP.* package is curated
by Bioconductor and maintains location and allele
information on SNP.

The hsSnpLocs environment is available through data(hsSnpLocs),
and contains a unified representation of the information in
the SNPlocs package (which is only available through various
chromosome-specific calls).
}
\value{two-row matrix -- top row is numeric suffix of
dbSNP ids, bottom row is location, genome-wide if chrnum is
missing
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(smlSet.example)
nn21 = colnames(smList(smlSet.example)[["21"]])
length(nn21)
ss = snpLocs.Hs(chrnum(21), rsid(nn21))
dim(ss)
# notes that the locations available are not completely
# congruent with those asserted in HapMap data
}
\keyword{ models }
