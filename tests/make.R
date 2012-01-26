library(GGBase)
if ("GGtools" %in% installed.packages()[,1]) {
 s20 = getSS("GGtools", "20")
 remk = make_smlSet( as(s20, "ExpressionSet"), smList(s20) )
 validObject(s20) & validObject(remk)
}
