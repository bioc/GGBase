
setClass("smlSet", contains="eSet",
   representation(smlEnv="environment", annotation="character",
     organism="character"),
   prototype=prototype(
       new("VersionedBiobase",
               versions=c(classVersion("eSet"), smlSet="1.1.3")),
           phenoData = new("AnnotatedDataFrame",
             data=data.frame(),
             varMetadata=data.frame(
               labelDescription=character(0))),
           annotation=character(0),
           organism=character(0),
           smlEnv = {e = new.env(); 
                     nl = list(chr1=matrix(), chr2=matrix()); 
                     assign("smList", nl, e); e},
           chromInds = numeric(0)))

setGeneric("smlEnv", function(x) standardGeneric("smlEnv"))
setMethod("smlEnv", "smlSet", function(x) x@smlEnv)
setGeneric("smList", function(x) standardGeneric("smList"))
setMethod("smList", "smlSet", function(x) x@smlEnv$smList)

# drop in Jan 2011
setMethod("exprs", "smlSet", function(object) {
  object@assayData$exprs}
)

valsml = function(object) {
 nexsamp = ncol(exprs(object))
 allns = sapply(smList(object), nrow)
 if (!(all(allns==allns[1])))
    return("varying numbers of rows in elements of smList")
 if (!(all(allns==nexsamp)))
    return("some SnpMatrix instances have nrows != ncol(exprs(smlSet))")
 if (is.null(names(smList(object))))
     return("smList elements must bear names e.g., c(1:22,'X', 'Y')")
 return(TRUE)
}

setValidity("smlSet", valsml)

setMethod("show", "smlSet", function(object) {
 cat("SnpMatrix-based genotype set:\n")
 cat("number of samples: ", length(sampleNames(object)), "\n")
 cat("number of chromosomes present: ", length(smList(object)), "\n")
 cat("annotation: ")
 cat( object@annotation, "\n" )
 if (length(dd <- dim(object@assayData$exprs))>0) {
  cat("Expression data dims:", dd[1], "x", dd[2], "\n")
 }
 nsnp = sum(sapply(smList(object), ncol))
  cat("Total number of SNP:", nsnp, "\n")
 cat("Phenodata: "); show(phenoData(object))
})


setMethod("[", "smlSet", function (x, i, j, ..., drop = FALSE) {
# j is strictly for samples
  if (!missing(j)) {
   # do snp matrices (samples are rows)
    L = smList(x)
# for snpMatrix2, omit drop spec
    LL = lapply(L, function(x) x[j,]) #,drop=FALSE] )
    ee = new.env()
    assign("smList", LL, ee)
    x@smlEnv = ee
   # do expression
    e = exprs(x)
    e = e[,j,drop=FALSE]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
   # do phenoData
    p = x@phenoData
    p = p[j,]
    x@phenoData = p
    protocolData(x) = protocolData(x)[j,]
  }
  if (!missing(i)) {
   if (is(i, "chrnum")) {
#
# odd use -- does not affect features -- but certainly does not
# affect samples -- it is an "assay" selection, so first coordinate seems ok
#
    L = smList(x)
    LL = L[i]
    ee = new.env()
    assign("smList", LL, ee)
    x@smlEnv = ee
    }
   else if (is(i, "probeId")) {
    e = exprs(x)
    e = e[i,,drop=FALSE]
    fd = featureData(x)[i,]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    x@featureData=fd
    }
   else if (is(i, "numeric")) {
    e = exprs(x)
    e = e[i,,drop=FALSE]
    fd = featureData(x)[i,]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    x@featureData=fd
    }
   else if (is(i, "genesym")) {
    e = exprs(x)
    e = e[sym2pid(x, i),,drop=FALSE]
    fd = featureData(x)[sym2pid(x, i),]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    x@featureData=fd
    }
   else if (is(i, "rsid")) {
    L = smList(x)
    nc = length(L)
    for (kk in 1:nc) {
      snin = colnames(L[[kk]])
      L[[kk]] = L[[kk]][, intersect(snin, i)]
      }
    ee = new.env()
    assign("smList", L, ee)
    x@smlEnv = ee
    }
   else stop(paste("[ method not defined for instance of ", class(i)))
  }
  return(x)
})

setMethod("combine", c("smlSet", "smlSet"), function(x, y, ...) {
 sx = smList(x)
 sy = smList(y)
 rsx = colnames(sx[[1]])
 rsy = colnames(sy[[1]])
 comm = intersect(rsx, rsy)
 sx = lapply(sx, function(x) x[, comm])
 sy = lapply(sy, function(x) x[, comm])
 fulls = list()
 for (i in 1:length(sx)) fulls[[i]] = rbind(sx[[i]], sy[[i]])
 names(fulls) = names(sx)
 ex = cbind(exprs(x), exprs(y))
 nad = assayDataNew("lockedEnvironment", exprs=ex)
 ee = new.env()
 assign("smList", fulls, ee)
 cx = x
 cx@assayData=nad
 cx@smlEnv = ee
 pdx = pData(x)
 pdy = pData(y)
 nnx = names(pdx)
 nny = names(pdy)
 cn = intersect(nnx, nny)
 pdat = rbind(pdx[,cn,drop=FALSE], pdy[,cn,drop=FALSE])
 pd = new("AnnotatedDataFrame", data=pdat)
 phenoData(cx) = pd
 cx
})


setAs("smlSet", "ExpressionSet", function(from) {
 ex = exprs(from)
 pd = phenoData(from)
 ans = new("ExpressionSet", exprs=ex, phenoData=pd)
 annotation(ans) = from@annotation
 ans
})

# private, important for plot_EvG

setGeneric("getAlleles", function(x, rs, ...) standardGeneric("getAlleles"))
setMethod("getAlleles", c("smlSet", "rsid"), function (x, rs)
{
    allrs = lapply(smList(x), colnames)
    #hits = lapply(allrs, function(x) grep(rs, x))
    hits = lapply(allrs, function(z) which(z == rs))
    kpi = sapply(hits, function(z) length(z) > 0)
    if (!(any(kpi))) stop("rs number not found in columns of smlSet")
    meta = list(chr = names(kpi[kpi]), col = hits[[which(kpi)]])
    ans = as(smList(x)[[meta$chr]][, meta$col], "character")
    if (any(ans == "Uncertain"))
    ans = as(smList(x)[[meta$chr]][, meta$col], "numeric")
    ans
})

