snpsNear = function (sym, radius = 1e+05, chrnum, ...) {
 if (is(sym, "GeneSet")) {
    if (geneIdType(sym)@type == "Annotation") alib = geneIdType(sym)@annotation
#    else if (geneIdType(sym)@type == "EntrezId") {
#      alib = "org.Hs.eg.db"
#      warning("assuming human organism for Entrez Id gene set")
#    }
    else stop("only Annotation-type Gene Sets handled at this time")
    sapply( geneIds(sym), function(x) snpsNear(probeId(x), radius=radius, annlib=alib))
    }
 else if (is(sym, "genesym") | is(sym, "probeId")) {
    pos = genePosition(sym, ...)
    chr = names(pos)
    if (chr %in% c("X", "Y")) 
        chr = ifelse(chr == "X", 23, 24)
    else chr = as.numeric(chr)
    ll = snpLocs.Hs(GGBase::chrnum(chr))
    inds = which(ll["loc",] >= pos-radius & ll["loc",] <= pos+radius)
    return(paste("rs", ll["rsid",inds], sep=""))
    }
 else if (is(sym, "rsid")) {
   allpos = snpLocs.Hs( chrnum(chrnum) )
   targ = allpos[ "loc",  allpos["rsid",] == as.numeric(gsub("rs", "", sym))]
   inds = which( allpos["loc",] >= targ - radius & allpos["loc",] <= targ + radius)
   return(paste("rs", allpos["rsid", inds], sep=""))
   }
 else if (is(sym, "numeric")) {
   targ = sym
   allpos = snpLocs.Hs( chrnum(chrnum) )
   inds = which( allpos["loc",] >= targ - radius & allpos["loc",] <= targ + radius)
   return(paste("rs", allpos["rsid", inds], sep=""))
   }
 else stop("need genesym or rsid or numeric instance")
}

setMethod("[", c("snp.matrix", "ANY", "rsid", "ANY"), 
  function (x, i, j, ..., drop = FALSE) 
  {
    cn = colnames(x)
    ii = intersect(cn, j)
    if (missing(i)) 
        x[, ii, drop = drop]
    else x[i, ii, drop = drop]
  })

