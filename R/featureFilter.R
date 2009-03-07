featureFilter = function(x, 
     requires=c("loc", "autosomal")) {
#
# will remove all features lacking keys in relevant CHRLOC, 
# having NA in CHRLOC (if "loc" %in% requires) , and having
# non autosomal address in CHRLOC (if "autosomal" %in% requires)
#
 clsuff = function(x) gsub(".db$", "", x)
 if (!is(x, "smlSet")) stop("requires smlSet input")
 fn = featureNames(x)
 require(paste(clsuff(Biobase::annotation(x)), ".db", sep=""), character.only=TRUE, quietly=TRUE)
 basic = keys(get(paste(clsuff(Biobase::annotation(x)), "CHRLOC", sep="")))
 bad = which(!(fn %in% basic))
 fn = fn[-bad]
 if ("loc" %in% requires) {
    locl = mget(fn, get(paste(clsuff(Biobase::annotation(x)), "CHRLOC", sep="")))
    locl = sapply(locl , "[", 1)
    bad = which(is.na(locl))
    fn = fn[-bad]
    }
 if ("autosomal" %in% requires){
    cn = mget(fn, get(paste(clsuff(Biobase::annotation(x)), "CHR", sep="")))
    cn = sapply(cn, "[", 1)
    bad = which(!(as.character(cn) %in% as.character(1:22)))
    fn = fn[-bad]
    }
 x[probeId(fn),]
}
