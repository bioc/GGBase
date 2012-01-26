library(GGBase)
if ("GGtools" %in% installed.packages()[,1]) {
 s20 = getSS("GGtools", "20")
 nsnp = ncol(smList(s20)[[1]])
 s20f = MAFfilter(s20, lower=.1)
 nsnp2 = ncol(smList(s20f)[[1]])
 ( nsnp == 119921 & nsnp2 == 46755 )
}
