
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

setMethod("plot", c("cwSnpScreenResult", "missing"),  # y bound to df request
  function(x, y=1, noSmooth=FALSE, npts=500, ...) {
   if (missing(y)) y = 1
   else if (!(y %in% c(1,2))) stop("df arg must be either 1 or 2")
   if (is(x@.Data[[1]], "snp.tests.glm"))  # for new approach, snpMatrix > 1.7
         allpv = p.value(x@.Data[[1]])
   else allpv = p.value(x@.Data[[1]]) #, y)
   kill = which(is.na(allpv))
   if (length(kill)>0) allpv = allpv[ -kill ]
   rsn = names(allpv)
   locstr = snpLocs.Hs(chrnum(x@chrnum), rsid(rsn)) # may not match all
   availRS = paste("rs", locstr["rsid",], sep="")
longnsubset = function (x, y) 
{
    mm = match(y, names(x))
    x[mm]
}

#   allpv = allpv[availRS]
   allpv = longnsubset(allpv, availRS)
   loc = locstr["loc",]
   if (noSmooth) plotf=plot
     else plotf=smoothScatter
   if (length(grep("resid", x@testType))>0) main = paste("resid", x@gene)
   else main=x@gene
   plotf(loc, -log10(allpv), main=main,
     xlab=paste("position on chr", x@chrnum),
     ylab=paste("-log10 p Gaussian LM [", y, "df]", sep=""), pch=19, cex=.8, ...)
   if (isCis(x)) 
        axis(3, at=genePosition(x@gene, annlib=x@annotation), col="red", lwd=2, label=" ")
})

#setMethod("plot", c("cwSnpScreenResult", "missing"),
#  function(x, y, noSmooth=FALSE, npts=500, ...) {
#   plot(x, y=1)
#})

#setMethod("plot", c("cwSnpScreenResult", "logical"),
#  function(x, y, noSmooth=FALSE, npts=500, ...) {
#   plot(x, y=1, noSmooth=noSmooth)
#})

setMethod("show", "multiGwSnpScreenResult", function(object) {
 cat("multi genome-wide snp screen result:\n")
 cat("gene set used as response:\n")
 show(object@geneset)
 cat("there are", length(object), "results.\n")
 cat("the call was:\n")
 print(object@call)
})
setMethod("show", "filteredMultiGwSnpScreenResult", function(object) {
 cat("filtered ")
 callNextMethod()
})
setMethod("show", "filteredGwSnpScreenResult", function(object) {
 cat("filtered ")
 callNextMethod()
})

setMethod("plot", c("cwSnpScreenResult", "character"),  # y bound to location package
  function(x, y="SNPlocs.Hsapiens.dbSNP.20090506", noSmooth=FALSE, npts=500, ...) {
   if (is(x@.Data[[1]], "snp.tests.glm"))  # for new approach, snpMatrix > 1.7
         allpv = p.value(x@.Data[[1]])
   else allpv = p.value(x@.Data[[1]]) #, y)
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

