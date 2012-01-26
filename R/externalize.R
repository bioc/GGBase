
externalize = function(smlSet, packname, author="Replace Me <auth@a.b.com>",
  maintainer="Replace Me <repl@a.b.com>") {
# creates folder structure for package
# saves expression data as ex in data folder
 dir.create(packname)
 dir.create(paste(packname, "/inst", sep=""))
 dir.create(paste(packname, "/R", sep=""))
 dir.create(paste(packname, "/inst/parts", sep=""))
 cn = names(smList(smlSet))
 partfol = paste(packname, "inst/parts/", sep="/")
 datfol = paste(packname, "data/", sep="/")
 for (i in cn) { assign(i, smList(smlSet)[[i]]); save(list=i, file=paste(partfol, i, ".rda", sep="")) }
 dir.create(paste(packname, "/data", sep=""))
 ex = as(smlSet, "ExpressionSet")
 save(ex, file=paste(datfol, "eset.rda", sep=""))
 dd = readLines(system.file("extpacksupp/DESCRIPTION.proto", package="GGtools"))
 zz = readLines(system.file("extpacksupp/zzz.R", package="GGtools"))
 dd = gsub("@MAINTAINER@", maintainer, dd)
 dd = gsub("@AUTHOR@", author, dd)
 dd = gsub("@PKGNAME@", packname, dd)
 writeLines(dd, paste(packname, "/DESCRIPTION", sep=""))
 writeLines(zz, paste(packname, "/R/zzz.R", sep=""))
 writeLines("", paste(packname, "/NAMESPACE", sep=""))
 cat(paste("now install", packname, "\n"))
 NULL
}

