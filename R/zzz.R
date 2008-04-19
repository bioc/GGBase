
.onLoad = function(libname, pkgname) {
 cat(paste(pkgname, "loading...\n"))
 require("ncdf", quietly=TRUE)
 ncpatt = ".*\\.nc$"
 allncdf = dir(system.file("extdata", package=pkgname), patt=ncpatt,
    full=TRUE)
 allsuff = dir(system.file("extdata", package=pkgname), patt=ncpatt)
 for (i in 1:length(allsuff)) {
   obname = paste(pkgname, allsuff[i], sep="_")
   AnnotationDbi:::addToNamespaceAndExport(obname, open.ncdf(allncdf[i]),
       pkgname)
   }
}

.onUnload <- function(libpath) {
 obs = objects("package:GGBase", patt=".*\\.nc$")
 sapply(obs, close.ncdf)
}

