\name{hmceuAmbB36_dbconn}
\alias{hmceuAmbB36_dbconn}
\alias{hmceuAmbB36_dbfile}
\title{SQLite infrastructure for SNP locations}
\description{SQLite infrastructure for SNP locations}
\usage{
hmceuAmbB36_dbconn()
hmceuAmbB36_dbfile()
}
\details{
hmceuAmbB36 refers to HapMap CEU Ambiguous IUPAC symbols
denoting SNP for hapmap build 36.

Follows conventions of SQLite annotation package
RDBMS referencing: the dbconn() function call above
evaluates to a SQLiteConnection instance; the dbfile() call
evaluates to a string giving the location of the SQLite file
holding the location table.
}
\value{
See details.
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
hmceuAmbB36_dbfile()
dbGetQuery(hmceuAmbB36_dbconn() , "select * from hmceuAmbB36 limit 5")
}
\keyword{ models }
