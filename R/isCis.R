
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
