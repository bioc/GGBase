snpsNear = function (gsym, radius = 1e+05) {
    pos = genePosition(gsym)
    chr = names(pos)
    if (chr %in% c("X", "Y")) 
        chr = ifelse(chr == "X", 23, 24)
    else chr = as.numeric(chr)
    ll = snpLocs.Hs(chrnum(chr))
    inds = which(ll["loc",] >= pos-radius & ll["loc",] <= pos+radius)
    paste("rs", ll["rsid",inds], sep="")
}

