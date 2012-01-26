.oldcode = function() {
#library(Biobase)
#else library(GSEABase)

# UPDATED AUGUST 2008 TO THIN OUT

# snpLocs -- use resources in Herve's SNPlocs package (with
# an interim step, see vignette newSnpLoc.Rnw

setClass("snpLocs", representation(locEnv="environment",
  offsets="numeric", organism="character", versions="character"))

#setMethod("show", "snpLocs", function(object) {
# cat("snpLocs instance, organism ", object@organism, "\n")
# cat("based on:\n")
# print(object@versions)
#})

# casting: chrnum

setClass("chrnum", contains="character")
setMethod("show", "chrnum", function(object) {
 cat("GGtools chrnum instance:\n")
 callNextMethod()
})
setGeneric("chrnum", function(x)standardGeneric("chrnum"))
setMethod("chrnum", "character", function(x)new("chrnum", x))
setMethod("chrnum", "numeric", function(x)new("chrnum", as.character(x)))

# casting: rsid

setClass("rsid", contains="character")
setMethod("show", "rsid", function(object) {
 cat("GGtools rsid instance:\n")
 callNextMethod()
})
setGeneric("rsid", function(x)standardGeneric("rsid"))
setMethod("rsid", "character", function(x)new("rsid", x))
setMethod("rsid", "numeric", function(x)new("rsid", as.character(x)))

# GGtools infrastructure from hm2ceu, Mar 31 2008 (c) VJ Carey
# revised Aug 2008


# bump class version to reflect accommodation of new eSet protocolData slot
setClass("smlSet", contains="eSet", 
   representation(smlEnv="environment", annotation="character",
     organism="character"),
   prototype=prototype(
       new("VersionedBiobase",
               versions=c(classVersion("eSet"), smlSet="1.1.2")),
           phenoData = new("AnnotatedDataFrame",
             data=data.frame(),
             varMetadata=data.frame(
               labelDescription=character(0))),
	   annotation=character(0),
	   organism=character(0),
#	   smlEnv = {e = new.env(); assign("smList", list(), e); e},  # this does not satisfy validity check
	   smlEnv = {e = new.env(); nl = list(chr1=matrix(), chr2=matrix()); assign("smList", nl, e); e},
           chromInds = numeric(0)))

setGeneric("smlEnv", function(x) standardGeneric("smlEnv"))
setMethod("smlEnv", "smlSet", function(x) x@smlEnv)
setGeneric("smList", function(x) standardGeneric("smList"))
setMethod("smList", "smlSet", function(x) x@smlEnv$smList)

# drop in Jan 2011
setMethod("exprs", "smlSet", function(object) {
  object@assayData$exprs}
)

valsml = function(object) {
 nexsamp = ncol(exprs(object))
 allns = sapply(smList(object), nrow)
 if (!(all(allns==allns[1]))) 
    return("varying numbers of rows in elements of smList")
 if (!(all(allns==nexsamp)))
    return("some SnpMatrix instances have nrows != ncol(exprs(smlSet))")
# if ((sl <- length(smList(object))) != (cl <- length(object@chromInds)))
#    return(paste("length of chromInds vector [", cl, "] not identical to that of smList(object) [",
#      sl, "]"))
# nna = names(annotation(object))
# if ((length(nna) != 2) || (!(all(nna == c("exprs", "snps")))))
#    return("annotation slot must be vector with names 'exprs' and 'snps'")
 if (is.null(names(smList(object))))
     return("smList elements must bear names e.g., c(1:22,'X', 'Y')")
 return(TRUE)
}

setValidity("smlSet", valsml)

# more casting

setClass("genesym", contains="character")
setGeneric("genesym", function(x) standardGeneric("genesym"))
setMethod("genesym", "character",  function(x) new("genesym", x))

setClass("probeId", contains="character")
setGeneric("probeId", function(x) standardGeneric("probeId"))
setMethod("probeId", "character",  function(x) new("probeId", x))

setClass("snpdepth", contains="numeric")
snpdepth = function(x) new("snpdepth", x)

setClass("phenoVar", contains="character")
setGeneric("phenoVar", function(x)standardGeneric("phenoVar"))
setMethod("phenoVar", "character", function(x)new("phenoVar", as.character(x)))

setClass("smlSummary", contains="list")
setMethod("show", "smlSummary", function(object) {
 cat("smList summary of length", length(object), "\n")
 cat("excerpt:\n")
 print(object[[1]][1:3,])
 cat("\n")
})

setClass("multiCisTestResult", representation(conditions="list",
   call="call"), contains="list")
setMethod("show", "multiCisTestResult", function(object) {
 cat("GGBase multiCisTestResult container.\n")
 cat("Call was:\n")
 print(object@call)
 cat("There are ", length(object), "results.\n")
 cat("Conditions raised for", sum(sapply(object@conditions, function(x)
     !is.na(x$cond))), "genes.\n")
})

setMethod("updateObject", "smlSet", function(object, ..., verbose=FALSE) {
  make_smlSet(as(object, "ExpressionSet"), smList(object), ...)
})

#setClass("fsmlSet", contains="smlSet")
#
#setAs("smlSet", "fsmlSet", function(from, to) {
#  z = smList(from)
#  fns = paste(substitute(from), "_", names(z), ".ff", sep="")
#  ffreflist = lapply( 1:length(z), function(w)
#     ff(initdata=z[[w]]@.Data, dim=dim(z[[w]]), vmode="raw", overwrite=TRUE,
#          dimnames=NULL, filename=fns[w]))
#  names(ffreflist) = names(z)
#  smlEnv = {e = new.env(); assign("smList", ffreflist, e); e}
#  new(to, smlEnv=smlEnv, annotation=from@annotation, organism=from@organism,
#        assayData=from@assayData, phenoData=from@phenoData, featureData=from@featureData,
#        experimentData=from@experimentData, protocolData=from@protocolData)
#})

#setClass("rangedSmlSet", representation(snplocs="GRangesList",
#    probelocs="GRangesList"), contains="smlSet")
#
##setMethod("initialize", "rangedSmlSet", function(.Object, ...) {
## .Object <- callNextMethod()
## .Object@snplocs = GRangesList()
## .Object@probelocs = GRangesList()
## .Object
##})
#
#setMethod("show", "rangedSmlSet", function(object) {
# callNextMethod()
# cat("rangeLists for SNP and probe locations available.\n")
#})
# 
#setGeneric("probeLocs<-", function(object, value)standardGeneric("probeLocs<-"))
#setMethod("probeLocs<-", c("smlSet", "GRangesList"), function(object, value) {
#   new("rangedSmlSet", snplocs=GRangesList(), probelocs=GRangesList(),
#       object) 
#})
 

genePosition = function (tok, genomeWide=FALSE, eglib="org.Hs.eg.db", annlib=NULL)
{
    tst = try(library(eglib, character.only=TRUE))
    egpref = gsub(".db", "", eglib)
    egSYM = get(paste(egpref, "SYMBOL", sep=""))
    egCHR = get(paste(egpref, "CHR", sep=""))
    egCHRLOC = get(paste(egpref, "CHRLOC", sep=""))
    egCHRLENGTHS = get(paste(egpref, "CHRLENGTHS", sep=""))
    if (!inherits(tst, "try-error")) {
        if (is(tok, "genesym")) 
            rmap = revmap(egSYM)
        else if (is(tok, "probeId")) {
            if (is.null(annlib)) stop("you must supply an annotation .db library name in annlib for a probeId token")
            require(annlib, character.only = TRUE, 
                quietly = TRUE)
            rmap = get(paste(gsub(".db", "", annlib), 
                "ENTREZID", sep = ""))
        }
        else {
            warning("tok is neither symbol nor probeID, we do not plot the location.")
            return(invisible(NULL))
        }
        egid = get(tok, rmap)
        ch = get(egid, egCHR)
        loc = get(egid, egCHRLOC)
        if (length(loc) > 1) {
            loc = loc[as.character(ch)][1]
        }
        if (length(loc) != 1) {
            warning("org.*.egCHRLOC has uninterpretable information for this gene; no tick at top of plot attempted.")
            return(invisible(NULL))
        }
        if (!genomeWide) return(abs(loc))
        chrl = egCHRLENGTHS
        chrbnd = cumsum(c(0,as.double(chrl[-length(chrl)])))
        if (ch == "X") 
            ch = 23
        else if (ch == "Y") 
            ch = 24
        else if (ch %in% c("Un", "MT")) 
            ch = 25
        else ch = as.numeric(ch)
        gpos = chrbnd[ch] + abs(loc)
        return(gpos)
    }
    warning("need org.*.eg.db for gene position.  no tick at top of plot attempted")
    return(invisible(NULL))
}

setMethod("plot", c("cwSnpScreenResult", "missing"),  # force y bound to location package
   function(x, y, noSmooth=FALSE, npts=500, ...) plot(x, "SNPlocs.Hsapiens.dbSNP.20090506",
     noSmooth, npts, ...))

setMethod("plot", c("cwSnpScreenResult", "character"),  # y bound to location package
  function(x, y="SNPlocs.Hsapiens.dbSNP.20090506", noSmooth=FALSE, npts=500, ...) {
   allpv = p.value(x@.Data[[1]])
   rsn = x@.Data[[1]]@snp.names
   names(allpv) = rsn
   #if (is(x@.Data[[1]], "snp.tests.glm"))  # for new approach, snpMatrix2 > 1.1
   #      allpv = p.value(x@.Data[[1]])
   #else allpv = p.value(x@.Data[[1]]) #, y)
   kill = which(is.na(allpv))
   if (length(kill)>0) allpv = allpv[ -kill ]
   rsn = names(allpv)
   SNY = do.call(":::", list(y, "SEQNAMES"))
   if (!(x@chrnum %in% SNY))  x@chrnum = chrnum(paste("chr", x@chrnum, sep=""))
   if (!(x@chrnum %in% SNY))  x@chrnum = chrnum(gsub("chr", "ch", x@chrnum))
   if (!(x@chrnum %in% SNY))  stop("attempts to harmonize @chrnum of cwSnpScreenResult object with SEQNAMES of y failed")
   loc = snpLocs.Hsapiens(rsn, x@chrnum, y) # may not match all
#   availRS = paste("rs", locstr["rsid",], sep="")
longnsubset = function (x, y)
{
    mm = match(y, names(x))
    x[mm]
}

   allpv = allpv[intersect(names(loc), names(allpv))]
   loc = loc[intersect(names(loc), names(allpv))]
#   allpv = longnsubset(allpv, availRS)
#   loc = locstr["loc",]
   if (noSmooth) plotf=plot
     else plotf=smoothScatter
   if (length(grep("resid", x@testType))>0) main = paste("resid", x@gene)
   else main=x@gene
   plotf(loc, -log10(allpv), main=main,
     xlab=paste("position on chr", x@chrnum),
     ylab=paste("-log10 p Gaussian LM [", 1, "df]", sep=""), pch=19, cex=.8, ...)
   if (isCis(x))
        axis(3, at=genePosition(x@gene, annlib=x@annotation), col="red", lwd=2, label=" ")
})


setGeneric("cwPlot", function(x,y,addloc,noSmooth,npts,...) standardGeneric("cwPlot"))
setMethod("cwPlot", c("cwSnpScreenResult", "character", "numeric"),  # y bound to location package
  function(x, y="SNPlocs.Hsapiens.dbSNP.20090506", addloc=NULL, noSmooth=FALSE, npts=500, ...) {
   if (is(x@.Data[[1]], "snp.tests.glm"))  # for new approach, snpMatrix > 1.7
         allpv = p.value(x@.Data[[1]])
   else allpv = p.value(x@.Data[[1]]) #, y)
   rsn = x@.Data[[1]]@snp.names
   names(allpv)= rsn
   kill = which(is.na(allpv))
   if (length(kill)>0) allpv = allpv[ -kill ]
   rsn = names(allpv)
   loc = snpLocs.Hsapiens(rsn, x@chrnum, y) # may not match all
#   availRS = paste("rs", locstr["rsid",], sep="")
longnsubset = function (x, y)
{
    mm = match(y, names(x))
    x[mm]
}

   if (!is.null(addloc)) loc = c(loc, addloc)
   allpv = allpv[intersect(names(loc), names(allpv))]
   loc = loc[intersect(names(loc), names(allpv))]
#   allpv = longnsubset(allpv, availRS)
#   loc = locstr["loc",]
   if (noSmooth) plotf=plot
     else plotf=smoothScatter
   if (length(grep("resid", x@testType))>0) main = paste("resid", x@gene)
   else main=x@gene
   plotf(loc, -log10(allpv), main=main,
     xlab=paste("position on chr", x@chrnum),
     ylab=paste("-log10 p Gaussian LM [", 1, "df]", sep=""), pch=19, cex=.8, ...)
   if (isCis(x))
        axis(3, at=genePosition(x@gene, annlib=x@annotation), col="red", lwd=2, label=" ")
})

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


#Class "cwSnpScreenResult"
#
#Slots:
#                                                                              
#Name:        .Data      chrnum        gene        psid  annotation    testType
#Class:        list      chrnum   character   character   character   character
#                              
#Name:         call sessionInfo
#Class:        call SessionInfo
#
isCis = function(x) {
 if (!is(x, "cwSnpScreenResult")) stop("only applies to cwSnpScreenResult")
 require(x@annotation, character.only=TRUE, quietly=TRUE)
 anm = gsub(".db$", "", x@annotation)
 ploc = get(x@psid, get(paste(anm, "CHR", sep="")))
 (ploc == as.character(x@chrnum) || paste("chr", x@chrnum, sep="") == ploc ||
    gsub("chr", "", x@chrnum) == ploc )
}

make_smlSet = function(es, sml, organism="Homo sapiens", harmonizeSamples=FALSE) {
 if (!inherits(es, "ExpressionSet")) stop("es must be ExpressionSet instance")
 if (!inherits(sml[[1]], "SnpMatrix")) stop("sml must be list of SnpMatrix instances from 'snpStats' package")
 if (is.null(names(sml))) stop("sml must be named list [typically list elements are named '1', '2', ... enumerating chromosomes, could just be 'all'")
 if (harmonizeSamples) {
  esn = sampleNames(es)
  ssn = rownames(sml[[1]]) # assumed common along list
  if (!all(esn %in% ssn) || !all(ssn %in% esn)) {
   warning("harmonizeSamples TRUE and sampleNames for es not coincident with rownames(sml[[1]]); harmonizing...")
   sn = intersect(esn, ssn)
   es = es[,sn]
   ns = names(sml)
   sml = lapply(sml, function(x)x[sn,])
   names(sml) = ns
   } 
 }  # otherwise validObject will get it
 smlenv = new.env()
 assign("smList", sml, envir=smlenv)
 new("smlSet", smlEnv=smlenv, annotation=annotation(es), organism=organism,
    assayData=assayData(es), phenoData=phenoData(es), 
    featureData=featureData(es), experimentData=experimentData(es))
}
 
pedinf2df = function(fn, ...) {
# create a data frame based on a hapmap pedigree info file
 V1 = V2 = V3 = V4 = V5 = V6 = V7 = 0
 x = read.table(fn, h=FALSE)
 attach(x)
 on.exit(detach(x))
 famid = V1
 persid = V2
 fathid = V3
 mothid = V4
 sampid = gsub("..*Sample:", "", as.character(V7))
 sampid = gsub(":.*$", "", sampid)
 isFounder = (V3 == 0 & V4 == 0)
 male = (V5 == 1)
 isAmom = persid %in% mothid
 isAdad = persid %in% fathid
 data.frame(famid=famid, persid=persid, mothid=mothid, fathid=fathid,
    sampid=sampid, isFounder=isFounder, male=male, isAmom = isAmom,
    isAdad = isAdad)
}
nsamp = function(x) ncol(exprs(x))

setMethod("show", "smlSet", function(object) {
 cat("SnpMatrix-based genotype set:\n")
 cat("number of samples: ", nsamp(object), "\n")
 cat("number of chromosomes present: ", length(smList(object)), "\n")
 cat("annotation: ")
 cat( object@annotation, "\n" )
 if (length(dd <- dim(object@assayData$exprs))>0) {
  cat("Expression data dims:", dd[1], "x", dd[2], "\n")
 }
 cat("Phenodata: "); show(phenoData(object))
})

#setMethod("exprs", "smlSet", function(object) {
#  object@assayData$exprs}
#)

setMethod("[", "smlSet", function (x, i, j, ..., drop = FALSE) {
# j is strictly for samples
  if (!missing(j)) {
   # do snp matrices (samples are rows)
    L = smList(x)
# for snpMatrix2, omit drop spec
    LL = lapply(L, function(x) x[j,]) #,drop=FALSE] )
    ee = new.env()
    assign("smList", LL, ee)
    x@smlEnv = ee
   # do expression
    e = exprs(x)
    e = e[,j,drop=FALSE]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
   # do phenoData
    p = x@phenoData
    p = p[j,]
    x@phenoData = p
    protocolData(x) = protocolData(x)[j,]
  }
  if (!missing(i)) {
   if (is(i, "chrnum")) {
#
# odd use -- does not affect features -- but certainly does not
# affect samples -- it is an "assay" selection, so first coordinate seems ok
#
    L = smList(x)
    LL = L[i]
    ee = new.env()
    assign("smList", LL, ee)
    x@smlEnv = ee
    }
   else if (is(i, "probeId")) {
    e = exprs(x)
    e = e[i,,drop=FALSE]
    fd = featureData(x)[i,]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    x@featureData=fd
    }
   else if (is(i, "numeric")) {
    e = exprs(x)
    e = e[i,,drop=FALSE]
    fd = featureData(x)[i,]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    x@featureData=fd
    }
   else if (is(i, "genesym")) {
    e = exprs(x)
    e = e[sym2pid(x, i),,drop=FALSE]
    fd = featureData(x)[sym2pid(x, i),]
    e = assayDataNew("lockedEnvironment", exprs=e)
    x@assayData = e
    x@featureData=fd
    }
   else if (is(i, "rsid")) {
    L = smList(x)
    nc = length(L)
    for (kk in 1:nc) {
      snin = colnames(L[[kk]])
      L[[kk]] = L[[kk]][, intersect(snin, i)]
      }
    ee = new.env()
    assign("smList", L, ee)
    x@smlEnv = ee
    }
   else stop(paste("[ method not defined for instance of ", class(i)))
  }
  return(x)
})


setGeneric("getAlleles", function(x, rs, ...) standardGeneric("getAlleles"))
setMethod("getAlleles", c("smlSet", "rsid"), function (x, rs)
{
    allrs = lapply(smList(x), colnames)
    #hits = lapply(allrs, function(x) grep(rs, x))
    hits = lapply(allrs, function(z) which(z == rs))
    kpi = sapply(hits, function(z) length(z) > 0)
    if (!(any(kpi))) stop("rs number not found in columns of smlSet")
    meta = list(chr = names(kpi[kpi]), col = hits[[which(kpi)]])
    ans = as(smList(x)[[meta$chr]][, meta$col], "character")
    if (any(ans == "Uncertain"))
    ans = as(smList(x)[[meta$chr]][, meta$col], "numeric")
    ans
})


setGeneric("plot_EvG", function(gsym, rsid, sms, ...) {
 standardGeneric("plot_EvG")})
setGeneric("plot_EvG2", function(gsym, rsid1, rsid2, sms, ...) {
 standardGeneric("plot_EvG2")})

sym2pid = function(sms, sym) {
  if (!is(sym, "genesym")) stop("sym2pid invoked without genesym instance")
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = get(sym, rmap)
  if (length(pid) > 1) {warning("multiple probes for this gene; taking first")}
  pid[1]
}

setMethod("plot_EvG", c("genesym", "rsid", "smlSet"),
 function(gsym, rsid, sms, ...) {
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = get(gsym, rmap)
  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
    an))
  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  thealleles = getAlleles(sms, rsid)
  gt = thealleles
  if (!is(thealleles[1], "numeric")) gt = factor(thealleles)
  if (is.factor(gt)) {
         plot(ex~gt, ylab=gsym, xlab=rsid, ...)
         points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
         } else {
         plot(ex~gt, ylab=gsym, xlab=paste("expected num. B alleles,", rsid), xlim=c(0,2), ...)
       }
  invisible(NULL)
})
setMethod("plot_EvG", c("probeId", "rsid", "smlSet"),
 function(gsym, rsid, sms, ...) {
#  an = sms@annotation
#  require(an, character.only=TRUE, quietly=TRUE)
#  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = gsym # get(gsym, rmap)
#  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
#    an))
#  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  thealleles = getAlleles(sms, rsid)
  gt = thealleles
  if (!is(thealleles[1], "numeric")) gt = factor(thealleles)
  if (is.factor(gt)) {
         plot(ex~gt, ylab=gsym, xlab=rsid, ...)
         points(jitter(as.numeric(gt),.4), ex, col="gray", pch=19)
         } else {
         plot(ex~gt, ylab=gsym, xlab=paste("expected num. B alleles,", rsid), xlim=c(0,2), ...)
       }
  invisible(NULL)
})

setMethod("plot_EvG2", c("genesym", "rsid", "rsid", "smlSet"),
 function(gsym, rsid1, rsid2, sms, ...) {
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = get(gsym, rmap)
  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
    an))
  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  gt1 = getAlleles(sms, rsid1)
  gt2 = getAlleles(sms, rsid2)
  if (is.character(gt1)) gt1 = factor(gt1)
  if (is.character(gt2)) gt2 = factor(gt2)
  if (is.numeric(gt1)) gtt = c(gt1,gt2)
  if (is.factor(gt1)) gtt = factor(paste(as.character(gt1), as.character(gt2)))
  plot(ex~gtt, ylab=gsym, xlab=paste("\n", rsid1, rsid2), ...)
  points(jitter(as.numeric(gtt),.4), ex, col="gray", pch=19)
})
setMethod("plot_EvG2", c("probeId", "rsid", "rsid", "smlSet"),
 function(gsym, rsid1, rsid2, sms, ...) {
#  an = sms@annotation
#  require(an, character.only=TRUE, quietly=TRUE)
#  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
#  pid = get(gsym, rmap)
  pid = gsym
#  if (length(pid) == 0) stop(paste("can't resolve", gsym, "in", 
#    an))
  if (length(pid) > 1) warning(paste("multiple probes for", gsym, "using first"))
  pid = pid[1]
  ex = exprs(sms)[pid,]
  gt1 = factor(getAlleles(sms, rsid1))
  gt2 = factor(getAlleles(sms, rsid2))
  gtt = factor(paste(as.character(gt1), as.character(gt2)))
  plot(ex~gtt, ylab=gsym, xlab=paste("\n", rsid1, rsid2), ...)
  points(jitter(as.numeric(gtt),.4), ex, col="gray", pch=19)
})

setGeneric("snpNames", function(x,c) standardGeneric("snpNames"))
setMethod("snpNames", c("smlSet", "missing"), function(x,c)
 lapply(smList(x), colnames))
setMethod("snpNames", c("smlSet", "chrnum"), function(x,c) {
 if (length(c)>1) stop("only scalar chrnum allowed")
 colnames(smList(x)[[c]])
})

setMethod("combine", c("smlSet", "smlSet"), function(x, y, ...) {
 sx = smList(x)
 sy = smList(y)
 rsx = colnames(sx[[1]])
 rsy = colnames(sy[[1]])
 comm = intersect(rsx, rsy)
 sx = lapply(sx, function(x) x[, comm])
 sy = lapply(sy, function(x) x[, comm])
 fulls = list()
 for (i in 1:length(sx)) fulls[[i]] = rbind(sx[[i]], sy[[i]])
 names(fulls) = names(sx)
 ex = cbind(exprs(x), exprs(y))
 nad = assayDataNew("lockedEnvironment", exprs=ex)
 ee = new.env()
 assign("smList", fulls, ee)
 cx = x
 cx@assayData=nad
 cx@smlEnv = ee
 pdx = pData(x)
 pdy = pData(y)
 nnx = names(pdx)
 nny = names(pdy)
 cn = intersect(nnx, nny)
 pdat = rbind(pdx[,cn,drop=FALSE], pdy[,cn,drop=FALSE])
 pd = new("AnnotatedDataFrame", data=pdat)
 phenoData(cx) = pd
 cx
})

setGeneric("snps", function(x, chr, ...) standardGeneric("snps"))
setMethod("snps", c("smlSet", "chrnum"), 
  function(x, chr, ...) {
     smList(x)[[chr]]
})
  

setAs("smlSet", "ExpressionSet", function(from) {
 ex = exprs(from)
 pd = phenoData(from)
 ans = new("ExpressionSet", exprs=ex, phenoData=pd)
 annotation(ans) = from@annotation
 ans
})


getSS = function( packname, chrs, renameChrs=NULL, probesToKeep=NULL,
   wrapperEndo=NULL ) {
 if (!is.null(renameChrs) && (length(chrs) != length(renameChrs)))
   stop("renameChrs must have same length as chrs in call to getSS")
 require(packname, character.only=TRUE)
 ex = get(load(system.file(package=packname, "data/eset.rda")))
 if (!is.null(probesToKeep)) ex = ex[probesToKeep,]
 partsfol = system.file("parts", package=packname)
 chk = sapply(chrs, function(x) file.exists(paste(partsfol, "/", x, ".rda", sep="")))
 if (!all(chk)) {
     cat("requesting ", paste(chrs, ".rda", sep=""))
     cat(" but finding\n")
     print(dir(partsfol))
     stop("cannot retrieve requested SNP file.")
     }
 sml = lapply(chrs, function(x) get(load(paste(partsfol, "/", x, ".rda", sep=""))))
 if (is.null(renameChrs)) names(sml) = chrs
 else names(sml) = renameChrs
 ans = make_smlSet( ex, sml, harmonizeSamples=TRUE )
 if (is.null(wrapperEndo)) return(ans) else {
   ans = wrapperEndo(ans)
   if (isTRUE(tst <- validObject(ans))) return(ans)
   stop(tst)
 }
}


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
#
#setGeneric("chromSpecLocs", function(sls, cnum, rsvec)
#  standardGeneric("chromSpecLocs"))
#
#setMethod("chromSpecLocs", c("snpLocs", "chrnum", "missing"),
#  function(sls, cnum, rsvec) {
#    get(cnum, sls@locEnv)
#})
#
#setMethod("chromSpecLocs", c("snpLocs", "chrnum", "character"),
#  function(sls, cnum, rsvec) {
#    rsvec = gsub("rs", "", rsvec)
#    cands = get(cnum, sls@locEnv)[1, ]
#    kp = match(as.integer(rsvec), cands, nomatch=0)
#    get(cnum, sls@locEnv)[, kp, drop=FALSE ]
#})
#
#setGeneric("genomeWideLocs", function(sls, rsvec)
#  standardGeneric("genomeWideLocs"))
#
#setMethod("genomeWideLocs", c("snpLocs", "character"),
#  function(sls, rsvec) {
#    rsvec = gsub("rs", "", rsvec)
#    allc = as.character(c(1:22, "X", "Y"))
#    allcand = NULL
#    off = sls@offsets
#    locs = rep(NA,length(rsvec))
#    rsid = rep(NA,length(rsvec))
#    chrn = rep(NA,length(rsvec))
#    cur = 1
#    for (i in 1:length(allc)) {
#       tmp = get(allc[i], sls@locEnv)
#       tmp2 = tmp[, match(rsvec, tmp[1,], nomatch=0), drop=FALSE]
#       nct= ncol(tmp2)
#       if (nct > 0) {
#           rsid[cur:(cur+nct-1)] = tmp2[1,]
#           locs[cur:(cur+nct-1)] = tmp2[2,]+off[i]
#           chrn[cur:(cur+nct-1)] = i
#           cur = cur+nct
#           }
#    }
#    rbind(rsid=rsid, gwloc=locs, chrn=chrn)
#})
#
#setGeneric("snpLocs.Hs", function(cnum, rsid) standardGeneric("snpLocs.Hs"))
#
#setMethod("snpLocs.Hs", c("chrnum", "rsid"), function(cnum, rsid) {
# if (!exists("hsSnpLocs")) data(hsSnpLocs)
# chromSpecLocs(hsSnpLocs, cnum, rsid)
#})
#
#setMethod("snpLocs.Hs", c("chrnum", "missing"), function(cnum, rsid) {
# if (!exists("hsSnpLocs")) data(hsSnpLocs)
# chromSpecLocs(hsSnpLocs, cnum)
#})
#
#setMethod("snpLocs.Hs", c("missing", "rsid"), function(cnum, rsid) {
# if (!exists("hsSnpLocs")) data(hsSnpLocs)
# genomeWideLocs(hsSnpLocs, rsid)
#})
# 
#setMethod("snpLocs.Hs", c("rsid"), function(cnum, rsid) {
# if (!exists("hsSnpLocs")) data(hsSnpLocs)
# genomeWideLocs(hsSnpLocs, cnum)
#})
#
#setGeneric("getSnpLocs", function(x,c) standardGeneric("getSnpLocs"))
#setMethod("getSnpLocs", c("smlSet", "missing"), function(x,c) {
# nn = snpNames(x)
# if (is(nn, "list")) {
#   ans = lapply(names(nn), function(x) snpLocs.Hs(chrnum(x)))
#   names(ans) = names(nn)
#   ans
#   }
# else snpLocs.Hs(rsid(nn))
#})
#
#setMethod("getSnpLocs", c("smlSet", "chrnum"), function(x,c) {
# nn = snpNames(x,c)
# return(snpLocs.Hs(rsid(nn)))
#})
snpLocs.Hsapiens = function( rsid, chrtok, spack="SNPlocs.Hsapiens.dbSNP.20090506" ) {
 require(spack, character.only=TRUE)
 ldf = getSNPlocs(chrtok)
 if (nrow(ldf) == 0) stop("chrtok must be wrong")
 locs = ldf$loc
 names(locs) = paste("rs", as.character(ldf$RefSNP_id), sep="")
 chk = intersect(rsid, names(locs))
 if (length(chk) != length(rsid)) warning(paste("some SNP in rsid were not found in location db", spack))
 locs[ chk ]
}
 
 
#snpsNear = function (sym, radius = 1e+05, chrnum, ...) {
# if (is(sym, "GeneSet")) {
#	stop("GeneSet inputs no longer supported")
##    if (geneIdType(sym)@type == "Annotation") alib = geneIdType(sym)@annotation
##    else if (geneIdType(sym)@type == "EntrezId") {
##      alib = "org.Hs.eg.db"
##      warning("assuming human organism for Entrez Id gene set")
##    }
##    else stop("only Annotation-type Gene Sets handled at this time")
##    return(sapply( geneIds(sym), function(x) try(snpsNear(probeId(x), radius=radius, annlib=alib))))
#    }
# else if (is(sym, "genesym") | is(sym, "probeId")) {
#    pos = genePosition(sym, ...)
#    what = as.numeric(pos)[1]
#    chr = names(pos)
#    if (chr %in% c("X", "Y")) 
#        chr = ifelse(chr == "X", 23, 24)
#    else chr = as.numeric(chr)
#    ll = snpLocs.Hs(GGBase::chrnum(chr))
#    inds = which(ll["loc",] >= pos-radius & ll["loc",] <= pos+radius)
#    ans = paste("rs", ll["rsid",inds], sep="")
#    attr(ans, "target") = c(chr=chr, loc=what)
#    return(ans)
#    }
# else if (is(sym, "rsid")) {
#   allpos = snpLocs.Hs( chrnum(chrnum) )
#   targ = allpos[ "loc",  allpos["rsid",] == as.numeric(gsub("rs", "", sym))]
#   what = targ
#   inds = which( allpos["loc",] >= targ - radius & allpos["loc",] <= targ + radius)
#   ans = paste("rs", allpos["rsid", inds], sep="")
#   attr(ans, "target") = c(chr=chrnum(chrnum), loc=what)
#   return(ans)
#   }
# else if (is(sym, "numeric")) {
#   targ = sym
#   allpos = snpLocs.Hs( chrnum(chrnum) )
#   inds = which( allpos["loc",] >= targ - radius & allpos["loc",] <= targ + radius)
#   ans = paste("rs", allpos["rsid", inds], sep="")
#   attr(ans, "target") = c(chr=chrnum(chrnum), loc=targ)
#   return(ans)
#   }
# else stop("need genesym or rsid or numeric instance")
#}

setMethod("[", c("SnpMatrix", "ANY", "rsid", "ANY"), 
  function (x, i, j, ..., drop = FALSE) 
  {
    cn = colnames(x)
    ii = intersect(cn, j)
    if (missing(i)) 
        x[, ii, drop = drop]
    else x[i, ii, drop = drop]
  })

# screening result containers

setOldClass("sessionInfo")
setClass("SessionInfo", contains="sessionInfo")
	   
setClass("gwSnpScreenResult", contains="list",
   representation(gene="character", psid="character", annotation="character", 
      testType="character", call="call",
        sessionInfo="SessionInfo")) #, modFmla="formula"))

setMethod("[", "gwSnpScreenResult", function(x, i,j,...,drop=FALSE) {
 if (!missing(j)) stop("no second subscript allowed")
 x@.Data = x@.Data[i]
 x@call = match.call()
 x
})

#setMethod("initialize", "gwSnpScreenResult", function(.Object, 
#     gene = character(0), psid = character(0), annotation=character(0),
#     testType=character(0), call=new("call"), modFmla=formula(a~b), ...) {
#     callNextMethod(.Object, gene=gene,
#		psid=psid, annotation=annotation, call=call, 
#		testType=testType, modFmla=modFmla, ...)
#})
#

setClass("cwSnpScreenResult", contains="gwSnpScreenResult",
   representation(chrnum="chrnum"))

#setMethod("initialize", "cwSnpScreenResult", function(.Object, ..., chrnum=new("chrnum")) {
#     .Object = callNextMethod(.Object, ...)
#     .Object@chrnum = chrnum
#     .Object
#})
# 

setMethod("show", "gwSnpScreenResult", function(object) {
 cat("gwSnpScreenResult for gene ", object@gene, " [probe ",
     object@psid, "]\n")
})

#setClass("multiGwSnpScreenResult", representation(geneset="GeneSet", call="call",
#   sessionInfo="SessionInfo"),
#   contains="list")
#setMethod("show", "cwSnpScreenResult", function(object) {
# cat("cwSnpScreenResult [chr", object@chrnum, "] for gene ", object@gene, " [probe ",
#     object@psid, "]\n")
#})

setClassUnion("cnumOrMissing", c("chrnum", "missing"))

#setClass("filteredGwSnpScreenResult", contains="gwSnpScreenResult")
#setClass("filteredMultiGwSnpScreenResult", contains="multiGwSnpScreenResult")


setMethod("[", "cwSnpScreenResult", function(x, i, j, ..., drop=FALSE) {
 d = x@.Data[[1]]  # must be length 1
 if (!missing(i)) {
   if (is(i, "numeric")) {
          d@chisq = d@chisq[i, , drop=FALSE]
          d@snp.names = d@snp.names[i]
          d@N = d@N[i]
          }
   else if (is(i, "character")) {
          kp = match(i, d@snp.names)
          kp = na.omit(kp)
	  d@snp.names = d@snp.names[kp]
          d@chisq = d@chisq[kp,,drop=FALSE]
	  d@N = d@N[kp]
          }
   }
 if (!missing(j)) stop("cannot select on 'column'")
 x@.Data[[1]] = d
 x
})

#setMethod("combine", c("multiGwSnpScreenResult",
#   "multiGwSnpScreenResult"), function(x, y, ...) {
#     d = c(x@.Data, y@.Data)
#     e = union(x@geneset, y@geneset)
#     m = match.call()
#     s = new("SessionInfo", sessionInfo())
#     new("multiGwSnpScreenResult", geneset=e, call=m,
#         sessionInfo=s, d)
#})
#
#setMethod("combine", c("filteredMultiGwSnpScreenResult",
#   "filteredMultiGwSnpScreenResult"), function(x, y, ...) {
#     d = c(x@.Data, y@.Data)
#     e = union(x@geneset, y@geneset)
#     m = match.call()
#     s = new("SessionInfo", sessionInfo())
#     new("filteredMultiGwSnpScreenResult", geneset=e, call=m,
#         sessionInfo=s, d)
#})


#.onLoad <- function(libname, package) {
#  library.dynam("GGBase", package)
#  methods:::bind_activation(TRUE)
#}
#
#.Last.lib <- function(libname, package) {
#  methods:::bind_activation(FALSE)
#}
}
