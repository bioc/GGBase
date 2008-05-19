\name{hmceuAmbB36_23a_dbconn}
\alias{hmceuAmbB36_23a_dbconn}
\alias{hmceuAmbB36_23a_dbfile}
\title{SQLite infrastructure for SNP locations}
\description{SQLite infrastructure for SNP locations}
\usage{
hmceuAmbB36_23a_dbconn()
hmceuAmbB36_23a_dbfile()
}
\details{
hmceuAmbB36 23a refers to HapMap CEU Ambiguous IUPAC symbols
denoting SNP for hapmap build 36 version 23a.

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
hmceuAmbB36_23a_dbfile()
con = hmceuAmbB36_23a_dbconn()
tab = dbListTables(con)
dbGetQuery(hmceuAmbB36_23a_dbconn() , paste("select * from", tab, " limit 5"))
}
\keyword{ models }
