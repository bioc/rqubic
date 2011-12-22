parseBic <- function(txt, delimiter=";") {
  features <- strsplit(txt[2], delimiter)[[1]]
  samples <- strsplit(txt[3], delimiter)[[1]]
  stopifnot(identical(txt[1],
                      paste(length(features), length(samples), sep=delimiter)))
  return(list(features=features, samples=samples))
}

readBiclusterResults <- function(filename,
                                 featureNames,
                                 sampleNames,
                                 delimiter=";",
                                 ...) {
  txt <- readLines(filename,...)
  stopifnot((length(txt)-1) %% 3L ==0)
  alg <- txt[1]
  bicno <- (length(txt)-1) %/% 3L
  bicfac <- factor(rep(seq(1L, bicno), each=3L))
  bil <- tapply(txt[2:length(txt)], bicfac, parseBic, delimiter=delimiter)
  
  uniqFeatures <- unique(sort(unlist(sapply(bil, function(x) x$features))))
  uniqSamples <- unique(sort(unlist(sapply(bil, function(x) x$samples))))
  
  if(!missing(featureNames)) {
    stopifnot(all(uniqFeatures %in% featureNames) && !anyDuplicated(featureNames))
    uniqFeatures <- featureNames
  }
  if(!missing(sampleNames)) {
    stopifnot(all(uniqSamples %in% sampleNames) && !anyDuplicated(sampleNames))
    uniqSamples <- sampleNames
  }
  
  nxc <- matrix(FALSE, nrow=bicno, ncol=length(uniqSamples))

  rxn <- sapply(1:bicno, function(i) uniqFeatures %in% bil[[i]]$features)
  rownames(rxn) <- uniqFeatures
  nxc <- t(sapply(1:bicno, function(i) uniqSamples %in% bil[[i]]$samples))
  colnames(nxc) <- uniqSamples
  
  obj <- new("Biclust",
             Parameters=list(Algorithm=alg),
             RowxNumber=rxn,
             NumberxCol=nxc,
             Number=bicno,
             info=list(featureNames=uniqFeatures, sampleNames=uniqSamples))
  return(obj)
}
