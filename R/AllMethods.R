##======================================================================##
## Methods of rqubic package
##
## Author: Jitao David Zhang <jitao_david.zhang@roche.com>, March 15, 2011
##======================================================================##

## Generics

##setGeneric("BCs", function(object) standardGeneric("BCs"))
##setGeneric("BCs<-", function(object,value) standardGeneric("BCs<-"))
setGeneric("BCcount", function(object) standardGeneric("BCcount"))
setGeneric("BCcount<-", function(object,value) standardGeneric("BCcount<-"))
setGeneric("Svalue", function(object, index) standardGeneric("Svalue"))
setGeneric("info", function(object,key, ...) standardGeneric("info"))
setGeneric("features", function(object) standardGeneric("features"))
setGeneric("featureCount", function(object, index) standardGeneric("featureCount"))
setGeneric("BCfeatures", function(object, index) standardGeneric("BCfeatures"))
setGeneric("BCfeatureCount", function(object, index) standardGeneric("BCfeatureCount"))
setGeneric("conditions", function(object) standardGeneric("conditions"))
setGeneric("conditionCount", function(object) standardGeneric("conditionCount"))
setGeneric("BCconditions", function(object, index) standardGeneric("BCconditions"))
setGeneric("BCconditionCount", function(object, index) standardGeneric("BCconditionCount"))
setGeneric("RowxNumber", function(object) standardGeneric("RowxNumber"))
setGeneric("NumberxCol", function(object) standardGeneric("NumberxCol"))
setGeneric("RowxNumber<-", function(object,value) standardGeneric("RowxNumber<-"))
setGeneric("NumberxCol<-", function(object,value) standardGeneric("NumberxCol<-"))
setGeneric("parameter", function(object, index) standardGeneric("parameter"))
setGeneric("quantileDiscretize", function(x,...) standardGeneric("quantileDiscretize"))
setGeneric("generateSeeds", function(object, ...) standardGeneric("generateSeeds"))

##----------------------------------------##
## helper functions
##----------------------------------------##
hasInputFeatureName <- function(x) RQUBIC_INPUT_FEATURENAME_SYMBOl %in% colnames(x)
inputFeatureName <- function(x) {
  if(RQUBIC_INPUT_FEATURENAME_SYMBOl %in% colnames(fData(x)))
    return(fData(x)[, RQUBIC_INPUT_FEATURENAME_SYMBOl])
  featureNames(assayData(x))
}
inputSampleName <- function(x) {
  if(RQUBIC_INPUT_SAMPLENAME_SYMBOl %in% colnames(fData(x)))
    fData(x)[, RQUBIC_INPUT_SAMPLENAME_SYMBOl]
  sampleNames(assayData(x))
}

##----------------------------------------##
## methods for rqubicSeeds
##----------------------------------------##
print.rqubicSeeds <- function(x,...) {
  cat("rqubic sorted seeds (associating C pointer):\n")
  cat(sprintf("  Number of seeds: %d\n", x[[1]]))
  cat(sprintf("  Minimum sample number: %d\n", attr(x, "minimumCol")))
}
minimumCol <- function(x) attr(x, "minimumCol")

##----------------------------------------##
## methods for Biclust
##----------------------------------------##
length.Biclust <- function(x) x@Number
dim.Biclust <- function(x) c(dim(RowxNumber(x))[1],
                             dim(NumberxCol(x))[2])

##length.QUBICBiclusterSet <- function(x) x@Number
##dim.QUBICBiclusterSet <- function(x) c(dim(RowxNumber(x))[1],
##                                       dim(NumberxCol(x))[2])

setMethod("RowxNumber", "Biclust", function(object) object@RowxNumber)
setMethod("NumberxCol", "Biclust", function(object) object@NumberxCol)
setReplaceMethod("RowxNumber", c("Biclust", "matrix"),
                 function(object, value) {
                   object@RowxNumber <- value
                   return(object)
                 })
setReplaceMethod("NumberxCol", c("Biclust", "matrix"),
                 function(object, value)  {
                 object@NumberxCol <- value
                 return(object)
               })

setMethod("BCcount", c("Biclust"), function(object) object@Number)
setMethod("BCcount<-", c("Biclust", "numeric"), function(object, value) object@Number <- as.integer(value))
setMethod("info", c("Biclust", "missing"), function(object, key) object@info)
setMethod("info", c("Biclust", "ANY"), function(object, key) object@info[[key]])
setMethod("Svalue", c("QUBICBiclusterSet", "missing"), function(object, index) info(object, "Svalues"))
setMethod("Svalue", c("QUBICBiclusterSet", "ANY"), function(object, index) Svalue(object)[index])

## features
setMethod("features", c("Biclust"), function(object) info(object, "features"))
setMethod("featureCount", c("Biclust"), function(object) length(features(object)))
setMethod("BCfeatures", c("Biclust", "missing"), function(object, index) {
  apply(RowxNumber(object),2,which)
})
setMethod("BCfeatures", c("Biclust", "ANY"), function(object, index) {
  res <- apply(RowxNumber(object),2,which)[index]
  if(length(res)==1)
    res <- res[[1]]
  return(res)
})
setMethod("BCfeatureCount", c("Biclust", "missing"), function(object) {
  colSums(RowxNumber(object))
})
setMethod("BCfeatureCount", c("Biclust", "ANY"), function(object, index) {
  res <- colSums(RowxNumber(object))[index]
  if(length(res)==1)
    res <- res[[1]]
  return(res)
})

## conditions
setMethod("conditions", c("Biclust"), function(object) info(object, "conditions"))
setMethod("conditionCount", c("Biclust"), function(object) length(conditions(object)))
setMethod("BCconditions", c("Biclust","missing"), function(object) {
  apply(NumberxCol(object), 1, which)
})
setMethod("BCconditions", c("Biclust", "ANY"), function(object, index) {
  res <- apply(NumberxCol(object), 1, which)[index]
  if(length(res)==1)
    res <- res[[1]]
  return(res)
})
setMethod("BCconditionCount", c("Biclust", "missing"), function(object) {
  rowSums(NumberxCol(object))
})
setMethod("BCconditionCount", c("Biclust", "ANY"), function(object, index) {
  res <- rowSums(NumberxCol(object))[index]
  if(length(res)==1)
    res <- res[[1]]
  return(res)
})

## features and conditions (compatible with ExpressionSet)
setMethod("featureNames", "Biclust", function(object) features(object))
setMethod("sampleNames", "Biclust", function(object) conditions(object))

## parameters
setMethod("parameter", c("Biclust", "missing"), function(object) object@Parameters)
setMethod("parameter", c("Biclust", "character"), function(object, index) {
  paras <- parameter(object)
  ind.match <- pmatch(index, names(paras), nomatch=NA_integer_, duplicates.ok=FALSE)
  if(length(index)!=1 || is.na(ind.match))
    stop("'index' must be one of the following parameters: ", paste(names(paras), collapse=","))
  paras[[ind.match]]
})

## subset and relevel
setMethod("[", c("Biclust", "ANY", "missing", "missing"),
          function(x, i, j, drop=FALSE) {
            new.obj <- x
            RowxNumber(new.obj) <- RowxNumber(x)[,i,drop=FALSE]
            NumberxCol(new.obj) <- NumberxCol(x)[i,,drop=FALSE]
            new.obj@Number <- ncol(RowxNumber(new.obj))
            return(new.obj)
          })


## show method modified from the Biclust package
setMethod("show", "QUBICBiclusterSet",
          function(object) {
            cat("\nAn object of class", class(object), "\n\n")
                        cat("Used features:", featureCount(object), "\n")
            cat("Used conditions:", conditionCount(object), "\n")
            cat(sprintf("Parameters: k=%d, f=%g, c=%g, o=%d, q=%g, r=%d\n",
                        parameter(object, "k"),
                        parameter(object, "f"),
                        parameter(object, "c"),
                        parameter(object, "o"),
                        parameter(object, "q"),
                        parameter(object, "r")))
            
            cat("Call:", deparse(object@Parameters$Call, 0.75 * getOption("width")), 
                sep = "\n\t")
            n <- BCcount(object)
            n <- min(c(n, 5))
            if (n > 1) {
              cat("\nNumber of Biclusters found: ", object@Number, "\n")
              cat("\nFirst ", n, " Cluster sizes:\n")
              rowcolsizes <- rbind(colSums(object@RowxNumber[, 1:n]), 
                                   rowSums(object@NumberxCol[1:n, ]))
              rownames(rowcolsizes) <- c("Number of Rows:", "Number of Columns:")
              colnames(rowcolsizes) <- paste("BC", 1:n)
              print(rowcolsizes)
            }
            else {
              if (n == 1) 
                cat("\nThere was one cluster found with\n ", sum(object@RowxNumber[, 
                                                                                   1]), "Rows and ", sum(object@NumberxCol), "columns")
              if (n == 0) 
                cat("\nThere was no cluster found")
            }

            cat("\n")
          })

##----------------------------------------##
## methods for other object types
##----------------------------------------##
setMethod("Svalue", c("matrix", "missing"), function(object) as.integer(prod(dim(object))))
setMethod("Svalue", c("eSet", "missing"), function(object) as.integer(prod(dim(object))))
