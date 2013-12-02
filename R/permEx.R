

.permEx = function(sms) {
#
# permute columns of expression data in an smlSet, ensuring that
# new numerical data are assigned to old sample names
#
# this is necessitated by snpStats::snp.rhs.tests, which will take
# steps to match response to predictor using sample id attributes
#
 ex = exprs(sms)
 nsamp = ncol(ex)
 pinds = sample(1:nsamp, size=nsamp, replace=FALSE)
 ini = colnames(ex)
 pex = ex[,pinds,drop=FALSE]
 colnames(pex) = ini  # this ensures that sample names are in the
                      # original order
 sms@assayData = assayDataNew("lockedEnvironment", exprs=pex)
 sms
}

setGeneric("permEx", function(sms)standardGeneric("permEx"))
setMethod("permEx", "smlSet", function(sms) .permEx(sms))
