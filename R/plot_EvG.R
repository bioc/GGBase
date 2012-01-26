
setGeneric("plot_EvG", function(gsym, rsid, sms, ...) {
 standardGeneric("plot_EvG")})

setMethod("plot_EvG", c("genesym", "rsid", "smlSet"),
 function(gsym, rsid, sms, ...) {
  an = sms@annotation
  require(an, character.only=TRUE, quietly=TRUE)
  rmap = revmap(get(paste(gsub(".db$", "", an), "SYMBOL", sep="")))
  pid = AnnotationDbi::get(gsym, rmap)
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
