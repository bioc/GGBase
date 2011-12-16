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
 if (length(bad)>0) fn = fn[-bad]
 if ("loc" %in% requires) {
    locl = mget(fn, clocOBJ, ifnotfound=NA)
    locln = lapply(locl, names)
    for (i in 1:length(locln)) {  # kick out hap/qbl references
       locl[[i]] = locl[[i]][ which(nchar(locln[[i]]) < 3) ]
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

