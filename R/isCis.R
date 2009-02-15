
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
 ploc = lookUp(x@psid, anm, "CHR")[[1]]
 ploc == as.character(x@chrnum)
}
