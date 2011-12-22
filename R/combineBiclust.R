## combine multiple Biclust class objects into one
setGeneric("combineBiclusts", function(x,y,...) standardGeneric("combineBiclusts"))

combine <- function(x,y,...) {
  biclist <- c(x,y, list(...))
  stopifnot(all(sapply(biclist, function(x) is(x, "Biclust"))))
  rns <- sapply(biclist, featureCount)
  cns <- sapply(biclist, conditionCount)
  stopifnot(identical(length(unique(rns)),1L))
  stopifnot(identical(length(unique(cns)),1L))

  new.rxn <- do.call(cbind, lapply(biclist, RowxNumber))
  new.nxc <- do.call(rbind, lapply(biclist, NumberxCol))
  new.number <- ncol(new.rxn)

  res <- biclist[[1]]
  res@RowxNumber <- new.rxn
  res@NumberxCol <- new.nxc
  res@Number <- new.number
  
  return(res)
}

combineQubic <- function(x,y,...) {
  biclist <- c(x,y, list(...))
  res <- combine(x,y,...)
  stopifnot(all(sapply(biclist, function(x) is(x, "QUBICBiclusterSet"))))
  res@Parameters$k <- sapply(biclist, function(x) parameter(x)$k)
  res@Parameters$f <- sapply(biclist, function(x) parameter(x)$f)
  res@Parameters$c <- sapply(biclist, function(x) parameter(x)$c)
  res@Parameters$o <- sapply(biclist, function(x) parameter(x)$o)
  res@Parameters$q <- sapply(biclist, function(x) parameter(x)$q)
  res@Parameters$r <- sapply(biclist, function(x) parameter(x)$r)
  res@info$Svalues <- unlist(lapply(biclist, Svalue))
  return(res)
}

setMethod("combineBiclusts",
          c("Biclust", "Biclust"), function(x,y,...) {
            combine(x,y,...)
          })

setMethod("combineBiclusts",
          c("QUBICBiclusterSet", "QUBICBiclusterSet"),
          function(x,y,...) {
            combineQubic(x,y,...)
          })
setMethod("combineBiclusts",
          c("list", "missing"), function(x,y,...) {
            do.call("combineBiclusts", x)
          })
