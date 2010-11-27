
make_smlSet = function(es, sml, organism="Homo sapiens") {
 if (!inherits(es, "ExpressionSet")) stop("es must be ExpressionSet instance")
 if (!inherits(sml[[1]], "SnpMatrix")) stop("sml must be list of SnpMatrix instances from 'snpMatrix2' package")
 if (is.null(names(sml))) stop("sml must be named list [typically list elements are named '1', '2', ... enumerating chromosomes, could just be 'all'")
 smlenv = new.env()
 assign("smList", sml, envir=smlenv)
 new("smlSet", smlEnv=smlenv, annotation=annotation(es), organism=organism,
    assayData=assayData(es), phenoData=phenoData(es), 
    featureData=featureData(es), experimentData=experimentData(es))
}
 
