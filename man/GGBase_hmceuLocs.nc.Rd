\name{GGBase_hmceuLocs.nc}
\alias{GGBase_hmceuLocs.nc}
\docType{data}
\title{ netCDF store (ncdf package) for SNP location data
for hapmap phase 2 SNP genotypes, build 36}
\description{ netCDF store (ncdf package) for SNP location data
for hapmap phase 2 SNP genotypes, build 36}
\usage{GGBase_hmceuLocs.nc}
\format{
  The format is that of an \code{\link[ncdf]{open.ncdf}} output.
}
\details{
This is a dynamically exported symbol that refers to an
S3 class ncdf instance, from ncdf package.
}
\source{
These locations are obtained in the 'support' data frames
returned from the snpMatrix read.HapMap.data calls.  The
reads occurred in early March 2008 and are from version 23;
 version 23a is available and will be imported soon.
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
library(GGBase)
GGBase_hmceuLocs.nc
}
\keyword{datasets}
