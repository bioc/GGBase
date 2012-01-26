
# casting: chrnum

setClass("chrnum", contains="character")
setMethod("show", "chrnum", function(object) {
 cat("GGBase chrnum instance:\n")
 callNextMethod()
})
setGeneric("chrnum", function(x)standardGeneric("chrnum"))
setMethod("chrnum", "character", function(x)new("chrnum", x))
setMethod("chrnum", "numeric", function(x)new("chrnum", as.character(x)))

# casting: rsid

setClass("rsid", contains="character")
setMethod("show", "rsid", function(object) {
 cat("GGBase rsid instance:\n")
 callNextMethod()
})
setGeneric("rsid", function(x)standardGeneric("rsid"))
setMethod("rsid", "character", function(x)new("rsid", x))
setMethod("rsid", "numeric", function(x)new("rsid", as.character(x)))

# more casting

setClass("genesym", contains="character")
setGeneric("genesym", function(x) standardGeneric("genesym"))
setMethod("genesym", "character",  function(x) new("genesym", x))

setClass("probeId", contains="character")
setGeneric("probeId", function(x) standardGeneric("probeId"))
setMethod("probeId", "character",  function(x) new("probeId", x))

