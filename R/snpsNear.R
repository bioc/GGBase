snpsNear = function (sym, radius = 1e+05, chrnum) {
 if (is(sym, "genesym")) {
    pos = genePosition(sym)
    chr = names(pos)
    if (chr %in% c("X", "Y")) 
        chr = ifelse(chr == "X", 23, 24)
    else chr = as.numeric(chr)
    ll = snpLocs.Hs(chrnum(chr))
    inds = which(ll["loc",] >= pos-radius & ll["loc",] <= pos+radius)
    return(paste("rs", ll["rsid",inds], sep=""))
    }
 else if (is(sym, "rsid")) {
   allpos = snpLocs.Hs( chrnum(chrnum) )
   targ = allpos[ "loc",  allpos["rsid",] == as.numeric(gsub("rs", "", sym))]
   inds = which( allpos["loc",] >= targ - radius & allpos["loc",] <= targ + radius)
   return(paste("rs", allpos["rsid", inds], sep=""))
   }
 else stop("need genesym or rsid instance")
}

