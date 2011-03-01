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

