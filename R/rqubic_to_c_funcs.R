eSetDimName <- function(eset, input, type=c("feature", "sample")) {
  if(!inherits(eset, "eSet"))
    stop("'eset' must be an object inheriting the 'eSet' class, for example an 'ExpressionSet'")
    
  type <- match.arg(type,
                    c("feature", "sample"))
  if(type=="feature") {
    names <- featureNames(eset)
    data <- fData(eset)
    rightDim <- dim(eset)[1]
  } else {
    names <- sampleNames(eset)
    data <- pData(eset)
    rightDim <- dim(eset)[2]
  }
  data.colnames <- colnames(data)
  
  if(!missing(input)) {
    if(length(input)==1 && is.character(input) && rightDim != 1) {
      if(!input %in% data.colnames) {
        errMsg <- paste("A character string as 'input' must match a column in fData/pData",
                        ifelse(dim(data)[2]==0, "",
                               paste(":\n ", paste(data.colnames, collapse=","))))
        stop(errMsg)
      }
      input <- as.character(data[, input])
    } else if(length(input) == rightDim) {
      input <- as.character(input)
    } else {
      errMsg <- paste("'input' Error, it can be",
                      "(*) a column name of fData/pData, or",
                      "(*) a vector giving names of all features/samples.",
                      "Otherwise its length must equal to the length of features/samples.",sep="\n")
      stop(errMsg)
    }
  } else {
    input <- names
  }
  return(input)
}

## write ExpressionSet to a matrix table required by the QUBIC command line tool
writeQubicInputFile <- function(x, file="",
                                featureNames, sampleNames) {
  ## check parameter sanities
  if(!inherits(x, "eSet"))
    stop("'x' must be an object inheriting the 'eSet' class, for example an 'ExpressionSet'")
  stopifnot(length(file)==1)

  ## following parameters are NOT passed onto the write.table: append, sep, quote, eol
  sep <- "\t"
  quote <- FALSE
  eol <- "\n"
  na <- "" ## TODO. check whether QUBIC accepts missing value as NA
  dec <- "."
  
  if(file=="") {
    file <- stdout()
  } else if (is.character(file)) {
    file <- file(file, "w")
    on.exit(close(file))
  } else if (!isOpen(file, "w")) {
    open(file, "w")
    on.exit(close(file))
  }
  if(!inherits(file, "connection"))
    stop("'file' must be a character string or connection")
  
  featureNames <- eSetDimName(x, featureNames, "feature")
  sampleNames <- eSetDimName(x, sampleNames, "sample")

  expmat <- exprs(x)
  rownames(expmat) <- featureNames
  
  ## qubic ask for an "o" appending before the colnames, no quote
  qubic.colnames <- c("o", sampleNames)
  writeLines(paste(qubic.colnames, collapse=sep), file, sep=eol)

  write.table(expmat, file=file, append=TRUE, quote=quote, sep=sep,
              eol=eol, na=na, dec=dec, row.names=TRUE, col.names=FALSE)
  
}
