
smlSummary = function(x) {
 if (!(is(x, "smlSet"))) stop("works only for smlSet instances")
 ans = lapply(smList(x), snpStats::col.summary)
 new("smlSummary", ans)
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
