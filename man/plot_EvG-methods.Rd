\name{plot_EvG-methods}
\docType{methods}
\alias{plot_EvG}
\alias{plot_EvG2}
\alias{plot_EvG-methods}
\alias{plot_EvG,genesym,rsid,smlSet-method}
\alias{plot_EvG,probeId,rsid,smlSet-method}
\alias{plot_EvG2,genesym,rsid,rsid,smlSet-method}
\alias{plot_EvG2,probeId,rsid,rsid,smlSet-method}
\title{ formal method for visualizing expression distributions vs genotype }
\description{
boxplot expression vs genotype
}
\section{Methods}{
\describe{

\item{gsym = "genesym", rsid = "rsid", sms = "smlSet"}{ generates an
annotated boxplot }
\item{multisnp methods}{plot\_EvG2 allows specification of
a second SNP rsid and shows boxplots over the
cross-tabulation of the allele combinations}
}}
\examples{
library(GGtools)
data(hmceuB36.2021)
plot_EvG(genesym("CPNE1"), rsid("rs6060535"), hmceuB36.2021)
}
\keyword{methods}
