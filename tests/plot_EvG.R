library(GGBase)
if ("GGtools" %in% installed.packages()[,1]) {
 s20 = getSS("GGtools", "20")
 plot_EvG(genesym("CPNE1"), rsid("rs6060535"), s20)
}
