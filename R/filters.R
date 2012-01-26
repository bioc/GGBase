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

MAFfilter = function(x, lower=0, upper=1) {
 if (!(is(x, "smlSet"))) stop("works only for smlSet instances")
 if (lower <= 0 & upper >= 1) return(x)
 ss = smlSummary(x)
 mafs = lapply(ss, "[", "MAF")
 allrs = lapply(ss, rownames)
 sml = smList(x)
 for (i in 1:length(mafs))
  {
  curok = which(mafs[[i]] >= lower & mafs[[i]] <= upper)
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

GTFfilter = function(x, lower=0) {
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

clipPCs = function (smlSet, inds2drop, center=TRUE)
{
#
# returns smlSet with transformed expressions --
# the principal components in inds2drop are omitted through
# zeroing components of the diagonal component of SVD of t(exprs)
#
    if (!is(smlSet, "smlSet"))
        stop("requires smlSet instance")
    ex = t(exprs(smlSet))
    ex = scale(ex, center=center, scale = FALSE)
    ss = svd(ex)
    d = ss$d
    d[inds2drop] = 0
    recon = t(ss$u %*% diag(d) %*% t(ss$v))
    rownames(recon) = featureNames(smlSet)
    colnames(recon) = sampleNames(smlSet)
    ne = assayDataNew("lockedEnvironment", exprs = recon)
    smlSet@assayData = ne
    smlSet
}

regressOut = function(sms, rhs, ...) {
 mm = model.matrix(rhs, data=pData(sms))
 f = limma::lmFit(exprs(sms), mm, ...)
 r = exprs(sms) - (f$coef %*% t(f$design))
 sms@assayData = assayDataNew("lockedEnvironment", exprs=r)
 sms
}

