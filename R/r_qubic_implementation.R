## check value(s) within lower and upper boundaries
isBound <- function(x, left, right, include.left=TRUE, include.right=TRUE) {
  stopifnot(length(left) ==1 && length(right) == 1)
  
  if(include.left)
    is.left <- x >= left
  else
    is.left <- x > left

  if(include.right)
    is.right <- x <= right
  else
    is.right <- x < right

  return(is.left & is.right)
}


##----------------------------------------------------------------------##
## disretize gene expression matrix by recursive quantlization of outliers
## S4 methods for matrices and ExpressionSets (eSet) 
##----------------------------------------------------------------------##

setMethod("quantileDiscretize", "matrix", function(x, q=0.06, rank=1L) {
  ## check q is a numeric value between 0 and 0.5
  is.valid.q <- length(q)==1 && is.numeric(q) && isBound(q, 0, 0.5, include.left=FALSE, include.right=FALSE)
  if(!is.valid.q)
    stop("'q' (quantile) must be a value in the range of (0,0.5)")

  ## rank is a positive integer larger than 1
  is.valid.rank <- length(rank)==1 && is.numeric(rank) && isBound(rank, left=0, right=dim(x)[2]/2, include.left=FALSE)
  if(!is.valid.rank) {
    stop("'rank' must be a positive integer less than half the count of conditions.")
  }
  rank <- as.integer(rank)
  
  ## the C function requires a REAL matrix, coercions
  storage.mode(x) <- "numeric"
  discr.x <- .Call("discretize_matrix", exp=x, q=q, rank=rank)

  return(discr.x)
})

setMethod("quantileDiscretize", "eSet", function(x, q=0.06, rank=1L) {
  x.exp <- exprs(x)
  x.exp.dis <- quantileDiscretize(x.exp, q, rank)
  exprs(x) <- x.exp.dis
  protocolData(x) <- new("AnnotatedDataFrame",
                              data=data.frame(quantile=q, rank=rank),
                              varMetadata=data.frame(labelDescription=c(q="quantile value", r="rank")))
  return(x)
})

## extract discretization parameters
discrParam <- function(eset, parameter=c("quantile", "rank")) {
  parameter <- match.arg(parameter)
  proData <- pData(protocolData(eset))
  if(!parameter %in% colnames(proData))
    return(NA)
  proData[, parameter]
}

##----------------------------------------------------------------------##
## generate seeds and sort them by coincidence scores
##----------------------------------------------------------------------##

setMethod("generateSeeds", "matrix", function(object, minColWidth=2L) {
  if(storage.mode(object) != "integer") {
    warning("'generateSeeds' requires an matrix of integer. Coerced\n")
    storage.mode(object) <- "integer"
  }

  is.valid.mincol <- length(minColWidth)==1 && minColWidth > 1L
  if(!is.valid.mincol)
    stop("'minColWidth' must be a positive integer")
  minColWidth <- as.integer(minColWidth)
  
  return(.Call("generate_sorted_seeds", object, minColWidth))
})

setMethod("generateSeeds", "eSet", function(object, minColWidth=2L) {
  object.exp <- exprs(object)
  object.seeds <- generateSeeds(object.exp, minColWidth=minColWidth)
  return(object.seeds)
})

##----------------------------------------------------------------------##
## Bicluster the seeds
##----------------------------------------------------------------------##
qubicRankSymbols <- function(rank=3L) {
  if(rank<1L)
    stop("rank must be a positive integer")
  rank <- as.integer(rank)
  c(0L, -(1L:rank),rank:1L)
}

rqubicVersionLabel <- function(rqubic.name="rqubic") {
  sprintf("%s %s", rqubic.name, package.version(rqubic.name))
}

## TODO: this function needs to be optimized for speed
blocksByIndex <- function(seeds, eset, index,
                          report.no,
                          tolerance,
                          filter.proportion) { 
  qubic.version <- rqubicVersionLabel()

  ## extract parameters 
  qubic.para.k <- minimumCol(seeds)
  qubic.para.f <- filter.proportion
  qubic.para.c <- tolerance
  qubic.para.o <- report.no
  qubic.para.q <- discrParam(eset, "quantile")
  qubic.para.r <- discrParam(eset, "rank")
  
  ## build qubic.bicObjects
  total.features <- featureNames(eset)
  total.samples <- sampleNames(eset)

  parameters <- list(Call=match.call(),
                     Method="rqubic",
                     qubicVersion=qubic.version,
                     inputDataFile=as.character(NA),
                     k=qubic.para.k,
                     f=qubic.para.f,
                     c=qubic.para.c,
                     o=qubic.para.o,
                     q=qubic.para.q,
                     r=qubic.para.r)
  number <- length(index)
  RowxNumber <- matrix(FALSE, nrow=length(total.features), ncol=number)
  NumberxCol <- matrix(FALSE, ncol=length(total.samples), nrow=number)
  for(i in seq(along=index)) {
    RowxNumber[index[[i]][[1]], i] <- TRUE
    NumberxCol[i, index[[i]][[2]]] <- TRUE
  }

  Svalues <- sapply(index, function(x) prod(length(x[[1]]), length(x[[2]])))
  
  obj <- new("QUBICBiclusterSet",
             Parameters=parameters,
             RowxNumber=RowxNumber,
             NumberxCol=NumberxCol,
             Number=number,
             info=list(features=total.features,
               conditions=total.samples,
               Svalues=Svalues))
  return(obj)
}

quBiclusterIndex <- function(seeds, eset, report.no=100L,
                             tolerance=0.95, filter.proportion=1) {
    ## parameter sanity checks
  is.valid.seeds <- inherits(seeds, "rqubicSeeds") & as.integer(seeds) > 0
  if(!is.valid.seeds)
    stop("'seeds' must be an object of rqubicSeeds, with at least one seed")

  exp <- exprs(eset)
  if(storage.mode(exp) != "integer") {
    warning("'eset' requires an matrix of integer as the value of 'exprs'. Coerced\n")
    storage.mode(exp) <- "integer"
  }
  
  is.valid.report.no <- is.numeric(report.no) && report.no > 0L
  if(!is.valid.report.no)
    stop("'report.no' must be a positive integer")
  report.no <- as.integer(report.no)

  is.valid.tolerance <- is.numeric(tolerance) && isBound(tolerance, 0, 1)
  if(!is.valid.tolerance) {
    warning("'tolerance' must be a positive number between 0 and 1. Set to default 0.95")
    tolerance <- 0.95
  }

  is.valid.filter <- is.numeric(filter.proportion) && isBound(filter.proportion, 0, 1)
  if(!is.valid.filter) 
    stop("'filter.proportion' must be a positive number between 0 and 1")
  storage.mode(filter.proportion) <- "double"

  sigma <- length(unique(as.vector(exp)))
  rank <- floor(sigma/2)
  symbols <- qubicRankSymbols(rank)

  ## level the exps
  if(!all(exp %in% symbols)) {
      stop("The exprs matrix of the ExpressionSet object must consist of only following integers:", paste(sort(symbols), collapse=","),"\n",
           "It is likely that the ExpressionSet has not been discretized. Please call quantileDiscretize first!")
  }
  exp.level <- factor(as.vector(exp), levels=symbols)
  exp.leveled <- matrix(as.integer(exp.level)-1L, ## -1L since C array starts with 0
                        nrow=nrow(exp), ncol=ncol(exp), byrow=FALSE)
  if(any(is.na(exp.leveled)))
      stop("Should not happen. Please contact the developer\n")

  ## rcIndex: an index of rows and columns indicating the bicluster memberships
  rcIndex <- .Call("qubicluster",
                   seeds=seeds,
                   exprs=exp.leveled,
                   sigma_val=sigma,
                   symbols=symbols,
                   report_no=report.no,
                   tolerance_val=tolerance,
                   filter_proportion=filter.proportion);
  
  return(rcIndex)
}

quBicluster <- function(seeds, eset, report.no=100L,
                        tolerance=0.95, filter.proportion=1) {
  ## rcIndex: an index of rows and columns indicating the bicluster memberships
  rcIndex <- quBiclusterIndex(seeds=seeds,
                              eset=eset,
                              report.no=report.no,
                              tolerance=tolerance,
                              filter.proportion=filter.proportion)
  
  ## subset expression with rcIndex
  obj <- blocksByIndex(seeds, eset, rcIndex,
                       report.no, tolerance, filter.proportion)
}
