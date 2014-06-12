
make_smlSet.legacy = function(es, sml, organism="Homo sapiens", harmonizeSamples=FALSE) {
 if (!inherits(es, "ExpressionSet")) stop("es must be ExpressionSet instance")
 if (!inherits(sml[[1]], "SnpMatrix")) stop("sml must be list of SnpMatrix instances from 'snpStats' package")
 if (is.null(names(sml))) stop("sml must be named list [typically list elements are named '1', '2', ... enumerating chromosomes, could just be 'all'")
 if (harmonizeSamples) {
  esn = sampleNames(es)
  ssn = rownames(sml[[1]]) # assumed common along list
  if (!all(esn %in% ssn) || !all(ssn %in% esn)) {
   message("harmonizeSamples TRUE and sampleNames for es not coincident with rownames(sml[[1]]); harmonizing...[not a warning any more]")
   sn = intersect(esn, ssn)
   es = es[,sn]
   ns = names(sml)
   sml = lapply(sml, function(x)x[sn,])
   names(sml) = ns
   } 
 }  # otherwise validObject will get it
 smlenv = new.env()
 assign("smList", sml, envir=smlenv)
 new("smlSet", smlEnv=smlenv, annotation=annotation(es), organism=organism,
    assayData=assayData(es), phenoData=phenoData(es), 
    featureData=featureData(es), experimentData=experimentData(es))
}
 
make_smlSet = function (es, sml, organism = "Homo sapiens", harmonizeSamples = FALSE)
{
    if (!inherits(es, "ExpressionSet"))
        stop("es must be ExpressionSet instance")
    if (!inherits(sml[[1]], "SnpMatrix"))
        stop("sml must be list of SnpMatrix instances from 'snpStats' package")
    if (is.null(names(sml)))
        stop("sml must be named list [typically list elements are named '1', '2', ... enumerating chromosomes, could just be 'all'")
    if (harmonizeSamples) {
        esn = sampleNames(es)
        ssn = rownames(sml[[1]])
        if (!all(esn %in% ssn) || !all(ssn %in% esn)) {
            message("harmonizeSamples TRUE and sampleNames for es not coincident with rownames(sml[[1]]); harmonizing...[not a warning]")
            sn = intersect(esn, ssn)
            es = es[, sn]
            ns = names(sml)
            sml = lapply(sml, function(x) x[sn, ])
            names(sml) = ns
        }
    }
    smlenv = new.env()
    assign("smList", sml, envir = smlenv)
    ans = new("smlSet")
    ans@smlEnv = smlenv
    ans@annotation = annotation(es)
    ans@organism = organism
    ans@assayData = assayData(es)
    ans@phenoData = phenoData(es)
    ans@protocolData = protocolData(es)
    ans@featureData = featureData(es)
    ans@experimentData = experimentData(es)
#    new("smlSet", smlEnv = smlenv, annotation = annotation(es),
#        organism = organism, assayData = assayData(es), phenoData = phenoData(es),
#        featureData = featureData(es), experimentData = experimentData(es))
    ans
}

