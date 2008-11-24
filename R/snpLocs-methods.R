
setGeneric("chromSpecLocs", function(sls, cnum, rsvec)
  standardGeneric("chromSpecLocs"))

setMethod("chromSpecLocs", c("snpLocs", "chrnum", "missing"),
  function(sls, cnum, rsvec) {
    get(cnum, sls@locEnv)
})

setMethod("chromSpecLocs", c("snpLocs", "chrnum", "character"),
  function(sls, cnum, rsvec) {
    rsvec = gsub("rs", "", rsvec)
    cands = get(cnum, sls@locEnv)[1, ]
    kp = match(as.integer(rsvec), cands, nomatch=0)
    get(cnum, sls@locEnv)[, kp, drop=FALSE ]
})

setGeneric("genomeWideLocs", function(sls, rsvec)
  standardGeneric("genomeWideLocs"))

setMethod("genomeWideLocs", c("snpLocs", "character"),
  function(sls, rsvec) {
    rsvec = gsub("rs", "", rsvec)
    allc = as.character(c(1:22, "X", "Y"))
    allcand = NULL
    off = sls@offsets
    locs = rep(NA,length(rsvec))
    rsid = rep(NA,length(rsvec))
    cur = 1
    for (i in 1:length(allc)) {
       tmp = get(allc[i], sls@locEnv)
       tmp2 = tmp[, match(rsvec, tmp[1,], nomatch=0), drop=FALSE]
       nct= ncol(tmp2)
       if (nct > 0) {
           rsid[cur:(cur+nct-1)] = tmp2[1,]
           locs[cur:(cur+nct-1)] = tmp2[2,]+off[i]
           cur = cur+nct
           }
    }
    rbind(rsid=rsid, loc=locs)
})

setGeneric("snpLocs.Hs", function(cnum, rsid) standardGeneric("snpLocs.Hs"))

setMethod("snpLocs.Hs", c("chrnum", "rsid"), function(cnum, rsid) {
 if (!exists("hsSnpLocs")) data(hsSnpLocs)
 chromSpecLocs(hsSnpLocs, cnum, rsid)
})

setMethod("snpLocs.Hs", c("chrnum", "missing"), function(cnum, rsid) {
 if (!exists("hsSnpLocs")) data(hsSnpLocs)
 chromSpecLocs(hsSnpLocs, cnum)
})

setMethod("snpLocs.Hs", c("missing", "rsid"), function(cnum, rsid) {
 if (!exists("hsSnpLocs")) data(hsSnpLocs)
 genomeWideLocs(hsSnpLocs, rsid)
})
 
setMethod("snpLocs.Hs", c("rsid"), function(cnum, rsid) {
 if (!exists("hsSnpLocs")) data(hsSnpLocs)
 genomeWideLocs(hsSnpLocs, cnum)
})

setGeneric("getSnpLocs", function(x,c) standardGeneric("getSnpLocs"))
setMethod("getSnpLocs", c("smlSet", "missing"), function(x,c) {
 nn = snpNames(x)
 if (is(nn, "list")) {
   ans = lapply(names(nn), function(x) snpLocs.Hs(chrnum(x)))
   names(ans) = names(nn)
   ans
   }
 else snpLocs.Hs(rsid(nn))
})

setMethod("getSnpLocs", c("smlSet", "chrnum"), function(x,c) {
 nn = snpNames(x,c)
 return(snpLocs.Hs(rsid(nn)))
})
