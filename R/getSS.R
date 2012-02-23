
getSS = function( packname, chrs, renameChrs=NULL, probesToKeep=NULL,
   wrapperEndo=NULL ) {
 if (!is.null(renameChrs) && (length(chrs) != length(renameChrs)))
   stop("renameChrs must have same length as chrs in call to getSS")
 #require(packname, character.only=TRUE)
 if (is.na(match(packname, installed.packages()[,1]))) stop(
		paste(packname, "not installed."))
 ex = get(load(system.file(package=packname, "data/eset.rda")))
 if (!is.null(probesToKeep)) ex = ex[probesToKeep,]
 partsfol = system.file("parts", package=packname)
 chk = sapply(chrs, function(x) file.exists(paste(partsfol, "/", x, ".rda", sep="")))
 if (!all(chk)) {
     cat("requesting ", paste(chrs, ".rda", sep=""))
     cat(" but finding\n")
     print(dir(partsfol))
     stop("cannot retrieve requested SNP file.")
     }
 sml = lapply(chrs, function(x) get(load(paste(partsfol, "/", x, ".rda", sep=""))))
 if (is.null(renameChrs)) names(sml) = chrs
 else names(sml) = renameChrs
 ans = make_smlSet( ex, sml, harmonizeSamples=TRUE )
 if (is.null(wrapperEndo)) return(ans) else {
   ans = wrapperEndo(ans)
   if (isTRUE(tst <- validObject(ans))) return(ans)
   stop(tst)
 }
}

