datacache <- new.env(hash=TRUE, parent=emptyenv())

hmceuAmbB36_dbconn <- function() dbconn(datacache) # lexical scope?
hmceuAmbB36_dbfile <- function() dbfile(datacache)

.onLoad <- function(libname, pkgname)
{
    require("methods", quietly=TRUE)
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "hmceuAmbB36.sql", 
       package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)
    ## Create the AnnObj instances
    aae = AnnotationDbi:::addToNamespaceAndExport
    aae("hmceuAmbB36_dbconn", hmceuAmbB36_dbconn, pkgname)
    aae("hmceuAmbB36_dbfile", hmceuAmbB36_dbfile, pkgname)
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(hmceuAmbB36_dbconn())
}

