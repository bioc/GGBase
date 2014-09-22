setClass("smlSummary", contains="list")
setMethod("show", "smlSummary", function(object) {
 cat("smList summary of length", length(object), "\n")
 cat("excerpt:\n")
 print(object[[1]][1:3,])
 cat("\n")
})

smlSummary = function(x) {
 if (!(is(x, "smlSet"))) stop("works only for smlSet instances")
 ans = lapply(smList(x), snpStats::col.summary)
 new("smlSummary", ans)
}


dropMonomorphies = function(sms) {
 sl = smList(sms)
 summs = lapply(sl, col.summary)
 todrop = lapply(summs, function(x) which(x[,"RAF"]==1 | x[,"RAF"]==0))
 for (i in 1:length(todrop))
   if (length(todrop[[i]])>0) sl[[i]] = sl[[i]][,-todrop[[i]]]
 sms@smlEnv$smList = sl
 sms
}

MAFfilter.legacy = function(x, lower=0, upper=1) {
 if (!(is(x, "smlSet"))) stop("works only for smlSet instances")
 if (lower <= 0 & upper >= 1) return(x)
 ss = smlSummary(x)
 mafs = lapply(ss, "[", "MAF")
 allrs = lapply(ss, rownames)
 sml = smList(x)
 for (i in 1:length(mafs))
  {
  curok = which(mafs[[i]] > lower & mafs[[i]] <= upper)
  if (length(curok) == 0) stop("limits eliminate all SNP on a chromosome, cannot proceed")
  kprs = allrs[[i]][curok]
  if (!all(allrs[[i]] %in% kprs))
     sml[[i]] = sml[[i]][, curok]
  }
 ne = new.env()
 assign("smList", sml, ne)
 x@smlEnv = ne
 x
}

MAFfilter = function (x, lower = 0, upper = 1)
{
    if (!(is(x, "smlSet")))
        stop("works only for smlSet instances")
    if (lower <= 0 & upper >= 1)
        return(x)
    sml <- x@smlEnv$smList
    maf = snpStats::col.summary(sml[[1]])[,"MAF",drop=FALSE]
    allrs = rownames(maf)
    curok = which(maf > lower & maf <= upper)
    rm(maf)
    if (length(curok) == 0)
            stop("limits eliminate all SNP on a chromosome, cannot proceed")
    if (length(curok) != length(allrs))
            x@smlEnv$smList[[1]] = x@smlEnv$smList[[1]][, curok]
    rm(allrs)
#    ne = new.env()
#    assign("smList", sml, ne)
#    x@smlEnv = ne
    x
}


GTFfilter.legacy = function(x, lower=0) {
 if (!(is(x, "smlSet"))) stop("works only for smlSet instances")
 if (lower <= 0 ) return(x)
 ss = smlSummary(x)
 mingtfs = lapply(ss, function(x) apply(x, 1, function(z) min(z[c("P.AA", "P.AB", "P.BB")])))
 allrs = lapply(ss, rownames)
 sml = smList(x)
 for (i in 1:length(mingtfs))
  {
  curok = which(mingtfs[[i]] >= lower)
  if (length(curok) == 0) stop("limits eliminate all SNP on a chromosome, cannot proceed")
  kprs = allrs[[i]][curok]
  if (!all(allrs[[i]] %in% kprs))
     sml[[i]] = sml[[i]][, curok]
  }
 ne = new.env()
 assign("smList", sml, ne)
 x@smlEnv = ne
 x
}

GTFfilter = function (x, lower = 0)
{
    if (!(is(x, "smlSet")))
        stop("works only for smlSet instances")
    if (lower <= 0)
        return(x)
    mingtf = rowMin(data.matrix(snpStats::col.summary(x@smlEnv$smList[[1]])[, c("P.AA", "P.AB", "P.BB")]))

    curok = which(mingtf >= lower)
    if (length(curok) == 0)  
            stop("limits eliminate all SNP on a chromosome, cannot proceed")
    if (!(length(curok) == length(mingtf)))
            x@smlEnv$smList[[1]] = x@smlEnv$smList[[1]][, curok]
    x
}


setMethod("nsFilter", "smlSet",
function (eset, require.entrez = TRUE, require.GOBP = FALSE,
    require.GOCC = FALSE, require.GOMF = FALSE, require.CytoBand = FALSE,
    remove.dupEntrez = TRUE, var.func = IQR, var.cutoff = 0.5,
    var.filter = TRUE, filterByQuantile = TRUE, feature.exclude = "^AFFX",
    ...) {
 ex = as(eset, "ExpressionSet")
 origfn = featureNames(ex)
 tmp = nsFilter(ex, require.entrez, require.GOBP,
      require.GOCC, require.GOMF, require.CytoBand,
         remove.dupEntrez, var.func, var.cutoff, var.filter, filterByQuantile,
        feature.exclude, ...)
 ex = tmp$eset
 ex = ex[ intersect(origfn, featureNames(ex)), ] # new defense for require.entrez=FALSE
 ans = make_smlSet(ex, smList(eset))
 experimentData(ans)@other = c(experimentData(ans)@other,
   nsfiltinfo = tmp[-1])
 ans
})




.clipPCs = function (smlSet, inds2drop, center=TRUE)
{
#
# returns smlSet with transformed expressions --
# the principal components in inds2drop are omitted through
# zeroing components of the diagonal component of SVD of t(exprs)
#
    if (!is(smlSet, "smlSet"))
        stop("requires smlSet instance")
    if (0 %in% inds2drop) {
       message("0 in inds2drop, so refraining from clipping; smlSet unchanged")
       return(smlSet)
       }
    ex = t(exprs(smlSet))
    recon = reconstruct( ex, inds2drop, center )
    rownames(recon) = featureNames(smlSet)
    colnames(recon) = sampleNames(smlSet)
    ne = assayDataNew("lockedEnvironment", exprs = recon)
    smlSet@assayData = ne
    smlSet
}

regressOut = function(sms, rhs, ...) {
 if (is(sms, "smlSet")) {
    mm = model.matrix(rhs, data=pData(sms))
    ex = exprs(sms)
    }
 else if (is(sms, "SummarizedExperiment")) {
    mm = model.matrix(rhs, data=colData(sms))
    message("using assay() to extract 'expression' matrix from SummarizedExperiment")
    ex = assay(sms)
    }
 else stop("only works for ExpressionSet or SummarizedExperiment")
 f = limma::lmFit(ex, mm, ...)
 r = ex - (f$coef %*% t(f$design))
 if (is(sms, "smlSet"))
     sms@assayData = assayDataNew("lockedEnvironment", exprs=r)
 else if (is(sms, "SummarizedExperiment"))
     assay(sms) = r
 sms
}

dropDupSNPs = function(sms, use.digest=TRUE, ...) {
 #require(digest)
 insml = smList(sms)
 if (length(insml)>1) stop("sms must have only one chromosome")
 sm = insml[[1]]
 if (use.digest) {
   cat("digesting...")
   dd = apply(sm@.Data,2,digest)
   }
 else {
   cat("coercing...")
   cgt = as(sm, "character")
   cat("pasting...")
   dd = apply(cgt,2,function(x)paste(x, collapse=":"))
   }
 cat("done.\n")
 dup = duplicated(dd)
 if (!any(dup)) return(sms)
 sm = sm[,-which(dup)]
 insml[[1]] = sm
 ne = new.env()
 assign("smList", insml, ne)
 sms@smlEnv = ne
 sms
}
# 
# 
#
setGeneric("clipPCs", 
 function(x, inds2drop, center=TRUE) standardGeneric("clipPCs"))

 
setMethod("clipPCs", 
  c("smlSet", "numeric", "logical"), function(x, inds2drop, center=TRUE){
   .clipPCs(smlSet=x, inds2drop, center)
})

setMethod("clipPCs", 
  c("SummarizedExperiment", "numeric", "logical"), function(x, inds2drop, center=TRUE){
   .clipPCs.SE(se=x, inds2drop, center)
})

setMethod("clipPCs", 
  c("SummarizedExperiment", "numeric", "missing"), function(x, inds2drop, center=TRUE){
   .clipPCs.SE(se=x, inds2drop, TRUE)
})

setMethod("clipPCs", 
  c("smlSet", "numeric", "missing"), function(x, inds2drop, center=TRUE){
   .clipPCs(smlSet=x, inds2drop, TRUE)
  })


reconstruct = function(ex, inds2drop, center=TRUE) {
    ex = scale(ex, center=center, scale = FALSE)
    ss = svd(ex)
    d = ss$d
    d[inds2drop] = 0
    t(ss$u %*% diag(d) %*% t(ss$v))
}

.clipPCs.SE = function(se, inds2drop, center=TRUE) {
     assn = names(assays(se))
     message(paste("clipping PCs", 
          paste0(selectSome(inds2drop),collapse=","), "from", assn[1], collapse=""))
     ex = t(assays(se)[[1]])
     recon = reconstruct(ex, inds2drop, center)
     assays(se)[[1]] = recon
     exptData(se)$PCsClipped = inds2drop
     se
     }

