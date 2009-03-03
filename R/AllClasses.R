#library(Biobase)
#else library(GSEABase)

# UPDATED AUGUST 2008 TO THIN OUT

# snpLocs -- use resources in Herve's SNPlocs package (with
# an interim step, see vignette newSnpLoc.Rnw

setClass("snpLocs", representation(locEnv="environment",
  offsets="numeric", organism="character", versions="character"))

setMethod("show", "snpLocs", function(object) {
 cat("snpLocs instance, organism ", object@organism, "\n")
 cat("based on:\n")
 print(object@versions)
})

# casting: chrnum

setClass("chrnum", contains="character")
setMethod("show", "chrnum", function(object) {
 cat("GGtools chrnum instance:\n")
 callNextMethod()
})
setGeneric("chrnum", function(x)standardGeneric("chrnum"))
setMethod("chrnum", "character", function(x)new("chrnum", x))
setMethod("chrnum", "numeric", function(x)new("chrnum", as.character(x)))

# casting: rsid

setClass("rsid", contains="character")
setMethod("show", "rsid", function(object) {
 cat("GGtools rsid instance:\n")
 callNextMethod()
})
setGeneric("rsid", function(x)standardGeneric("rsid"))
setMethod("rsid", "character", function(x)new("rsid", x))
setMethod("rsid", "numeric", function(x)new("rsid", as.character(x)))

# GGtools infrastructure from hm2ceu, Mar 31 2008 (c) VJ Carey
# revised Aug 2008

valsml = function(object) {
 allns = sapply(smList(object), nrow)
 if (!(all(allns==allns[1]))) 
    return("varying numbers of rows in elements of smList")
# if ((sl <- length(smList(object))) != (cl <- length(object@chromInds)))
#    return(paste("length of chromInds vector [", cl, "] not identical to that of smList(object) [",
#      sl, "]"))
# nna = names(annotation(object))
# if ((length(nna) != 2) || (!(all(nna == c("exprs", "snps")))))
#    return("annotation slot must be vector with names 'exprs' and 'snps'")
 if (is.null(names(smList(object))))
     return("smList elements must bear names e.g., c(1:22,'X', 'Y')")
 return(TRUE)
}

setClass("smlSet", contains="eSet", 
   representation(smlEnv="environment", annotation="character",
     organism="character"),
   validity=valsml, prototype=prototype(
       new("VersionedBiobase",
               versions=c(classVersion("eSet"), smlSet="1.1.1")),
           phenoData = new("AnnotatedDataFrame",
             data=data.frame(),
             varMetadata=data.frame(
               labelDescription=character(0))),
	   annotation=character(0),
	   organism=character(0),
	   smlEnv = {e = new.env(); assign("smList", list(), e); e},
           chromInds = numeric(0)))

setGeneric("smlEnv", function(x) standardGeneric("smlEnv"))
setMethod("smlEnv", "smlSet", function(x) x@smlEnv)
setGeneric("smList", function(x) standardGeneric("smList"))
setMethod("smList", "smlSet", function(x) x@smlEnv$smList)

# more casting

setClass("genesym", contains="character")
setGeneric("genesym", function(x) standardGeneric("genesym"))
setMethod("genesym", "character",  function(x) new("genesym", x))

setClass("probeId", contains="character")
setGeneric("probeId", function(x) standardGeneric("probeId"))
setMethod("probeId", "character",  function(x) new("probeId", x))

setClass("snpdepth", contains="numeric")
snpdepth = function(x) new("snpdepth", x)

setClass("phenoVar", contains="character")
setGeneric("phenoVar", function(x)standardGeneric("phenoVar"))
setMethod("phenoVar", "character", function(x)new("phenoVar", as.character(x)))
