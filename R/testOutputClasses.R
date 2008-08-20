# screening result containers
	   
setClass("gwSnpScreenResult", contains="list",
   representation(gene="character", psid="character", annotation="character", testType="character",
      call="call"))

setMethod("show", "gwSnpScreenResult", function(object) {
 cat("gwSnpScreenResult for gene ", object@gene, " [probe ",
     object@psid, "]\n")
})

setClass("multiGwSnpScreenResult", representation(geneset="GeneSet", call="call"),
   contains="list")

setClass("cwSnpScreenResult", contains="gwSnpScreenResult",
   representation(chrnum="chrnum"))
setMethod("show", "cwSnpScreenResult", function(object) {
 cat("cwSnpScreenResult [chr", object@chrnum, "] for gene ", object@gene, " [probe ",
     object@psid, "]\n")
})

setClassUnion("cnumOrMissing", c("chrnum", "missing"))

setClass("filteredGwSnpScreenResult", contains="gwSnpScreenResult")
setClass("filteredMultiGwSnpScreenResult", contains="multiGwSnpScreenResult")

