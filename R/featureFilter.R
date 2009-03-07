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
 clocOBJ = get(paste(clsuff(Biobase::annotation(x)), "CHRLOC", sep=""))
 cOBJ = get(paste(clsuff(Biobase::annotation(x)), "CHR", sep=""))
 basic = keys(clocOBJ)
 bad = which(!(fn %in% basic))
 fn = fn[-bad]
 if ("loc" %in% requires) {
    locl = mget(fn, clocOBJ, ifnotfound=NA)
    locln = lapply(locl, names)
    for (i in 1:length(locln)) {  # kick out hap/qbl references
       locl[[i]] = locl[[i]][ nchar(locln[[i]]) < 3 ]
    }
    locl = sapply(locl , "[", 1)  # take first qualifying location
    bad = which(is.na(locl))
    fn = fn[-bad]
    }
 if ("autosomal" %in% requires){
    cn = mget(fn, cOBJ, ifnotfound=NA)
    cn = sapply(cn, "[", 1)
    bad = which(!(as.character(cn) %in% as.character(1:22)))
    fn = fn[-bad]
    }
 x[probeId(fn),]
}
