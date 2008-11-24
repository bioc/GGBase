nsamp = function(x) ncol(exprs(x))

setMethod("show", "smlSet", function(object) {
 cat("snp.matrix-based genotype set:\n")
 cat("number of samples: ", nsamp(object), "\n")
 cat("number of snp.matrix: ", length(smList(object)), "\n")
 cat("annotation: ")
 cat( object@annotation, "\n" )
 if (length(dd <- dim(object@assayData$exprs))>0) {
  cat("Expression data dims:", dd[1], "x", dd[2], "\n")
 }
 cat("Phenodata: "); show(phenoData(object))
})

setMethod("exprs", "smlSet", function(object) {
  object@assayData$exprs}
)

setMethod("[", "smlSet", function (x, i, j, ..., drop = FALSE) {
  if (!missing(j)) {
   # do snps
    L = smList(x)
    LL = lapply(L, function(x) x[j,,drop=FALSE] )
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
  }
  if (!missing(i)) {
   if (is(i, "chrnum")) {
    L = smList(x)
    LL = L[i]
    ee = new.env()
    assign("smList", LL, ee)
    x@smlEnv = ee
    }
   else if (is(i, "probeId")) {
    e = exprs(x)
    e = e[i,,drop=FALSE]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    }
   else if (is(i, "numeric")) {
    e = exprs(x)
    e = e[i,,drop=FALSE]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    }
   else if (is(i, "genesym")) {
    e = exprs(x)
    e = e[sym2pid(x, i),,drop=FALSE]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    }
   else stop(paste("[ method not defined for instance of ", class(i)))
  }
  return(x)
})


setGeneric("getAlleles", function(x, rs, ...) standardGeneric("getAlleles"))
setMethod("getAlleles", c("smlSet", "rsid"), function (x, rs)
{
    allrs = lapply(smList(x), colnames)
    #hits = lapply(allrs, function(x) grep(rs, x))
    hits = lapply(allrs, function(z) which(z == rs))
    kpi = sapply(hits, function(z) length(z) > 0)
    if (!(any(kpi))) stop("rs number not found in columns of smlSet")
    ans = list(chr = names(kpi[kpi]), col = hits[[which(kpi)]])
    as(smList(x)[[ans$chr]][, ans$col], "character")
})


setGeneric("plot_EvG", function(gsym, rsid, sms, ...) {
 standardGeneric("plot_EvG")})
setGeneric("plot_EvG2", function(gsym, rsid1, rsid2, sms, ...) {
 standardGeneric("plot_EvG2")})

sym2pid = function(sms, sym) {
  if (!is(sym, "genesym")) stop("sym2pid invoked without genesym instance")
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = get(sym, rmap)
  if (length(pid) > 1) {warning("multiple probes for this gene; taking first")}
  pid[1]
}

setMethod("plot_EvG", c("genesym", "rsid", "smlSet"),
 function(gsym, rsid, sms, ...) {
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = get(gsym, rmap)
  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
    an))
  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  gt = factor(getAlleles(sms, rsid))
  plot(ex~gt, ylab=gsym, xlab=rsid, ...)
  points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
})
setMethod("plot_EvG", c("probeId", "rsid", "smlSet"),
 function(gsym, rsid, sms, ...) {
#  an = sms@annotation
#  require(an, character.only=TRUE, quietly=TRUE)
#  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = gsym # get(gsym, rmap)
#  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
#    an))
#  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  gt = factor(getAlleles(sms, rsid))
  plot(ex~gt, ylab=gsym, xlab=rsid, ...)
  points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
})

setMethod("plot_EvG2", c("genesym", "rsid", "rsid", "smlSet"),
 function(gsym, rsid1, rsid2, sms, ...) {
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = get(gsym, rmap)
  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
    an))
  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  gt1 = factor(getAlleles(sms, rsid1))
  gt2 = factor(getAlleles(sms, rsid2))
  gtt = factor(paste(as.character(gt1), as.character(gt2)))
  plot(ex~gtt, ylab=gsym, xlab=paste("\n", rsid1, rsid2), ...)
  points(jitter(as.numeric(gtt),.4), ex, col="gray", pch=19)
})
setMethod("plot_EvG2", c("probeId", "rsid", "rsid", "smlSet"),
 function(gsym, rsid1, rsid2, sms, ...) {
#  an = sms@annotation
#  require(an, character.only=TRUE, quietly=TRUE)
#  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
#  pid = get(gsym, rmap)
  pid = gsym
#  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
#    an))
  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  gt1 = factor(getAlleles(sms, rsid1))
  gt2 = factor(getAlleles(sms, rsid2))
  gtt = factor(paste(as.character(gt1), as.character(gt2)))
  plot(ex~gtt, ylab=gsym, xlab=paste("\n", rsid1, rsid2), ...)
  points(jitter(as.numeric(gtt),.4), ex, col="gray", pch=19)
})

setGeneric("snpNames", function(x,c) standardGeneric("snpNames"))
setMethod("snpNames", c("smlSet", "missing"), function(x,c)
 lapply(smList(x), colnames))
setMethod("snpNames", c("smlSet", "chrnum"), function(x,c) {
 if (length(c)>1) stop("only scalar chrnum allowed")
 colnames(smList(x)[[c]])
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

setGeneric("snps", function(x, chr, ...) standardGeneric("snps"))
setMethod("snps", c("smlSet", "chrnum"), 
  function(x, chr, ...) {
     smList(x)[[chr]]
})
  
