##======================================================================##
## Parsers for QUBIC ANSI C program outputs
##
## Author: Jitao David Zhang <jitao_david.zhang@roche.com>, March 15, 2011
##======================================================================##

parseQubicRules <- function(filename) {
  rules <- scan(filename, what="character", sep="\n", quiet=TRUE)

  ## the pattern needs to be revised!
  rule.pattern <- "^row\\s([[:graph:]]*)\\s:low=((\\d|\\.|-)*),\\sup=((\\d|\\.|-)*);\\s(\\d*)\\sdown-regulated,(\\d*)\\sup-regulated"
  rules.notfollow.pattern <- !grepl(rule.pattern, rules)
  if(any(rules.notfollow.pattern)) {
    warning.lines <-  paste("[", which(rules.notfollow.pattern), "]:", rules[rules.notfollow.pattern], collapse="\n", sep="")
    stop(sprintf("The following %d lines do not seem to be QUBIC rules:\n%s",
                 sum(rules.notfollow.pattern),
                 warning.lines))
  }
    
  rules.parsed <- gsub(rule.pattern, "\\1\t\\2\t\\4\t\\6\t\\7", rules)
  rules.split <- strsplit(rules.parsed, "\t")

  rules.mat <- do.call(rbind, rules.split)
  rules.df <- data.frame(probe=rules.mat[,1],
                         low=as.numeric(rules.mat[,2]),
                         up=as.numeric(rules.mat[,3]),
                         down.count=as.integer(rules.mat[,4]),
                         up.count=as.integer(rules.mat[,5]))
  return(rules.df)
}


parseQubicChars <- function(file, check.names=FALSE,...) {
  Call <- match.call(expand.dots=TRUE)
  for(argname in c("sep", "head", "row.names"))
    if(!is.null(Call[[argname]]))
      warning(gettextf("attempt to set '%s' ignored", argname), domain=NA)
  Call$sep <- "\t"
  Call$head <- TRUE
  Call$check.names <- check.names
  Call[[1L]] <- as.name("read.csv")
  res <- eval(Call)
  ## note that the first column cannot be set to row names since they might duplicate
  colnames(res)[1] <- "probe"
  return(res)
}

parseQubicCountedVector <- function(vector) {
  bc.items <- strsplit(vector, " ")[[1]]
  bc.item.count <- as.integer(gsub("\\[(\\d*)\\]:", "\\1", bc.items[[3]]))
  bc.items <- bc.items[4:length(bc.items)]
  if(is.na(bc.item.count) || bc.item.count != length(bc.items))
    stop("The count in the vector does not equal the count reported. File may be corrupted")

  return(bc.items)
}


## TODO: the QubicBiclusterSet object should use index to mark biclusters, but not eSets: too heavy
parseQubicBicluster <- function(vector) {
  bc.head.pattern <- "^BC\\d*\\sS=(\\d*)$"
  bc.Svalue <- as.integer(gsub(bc.head.pattern, "\\1", vector[1]))

  ## bicluster features
  bc.features <- parseQubicCountedVector(vector[[2]])

  ## bicluster samples
  bc.conditions <- parseQubicCountedVector(vector[[3]])

  ## table
  vector.split <- strsplit(vector[-(1:3)], "\t")
  ## check features/cond correspond
  table.features.raw <- sapply(vector.split, "[[", 1)
  table.features <- gsub("\\s*([[:graph:]]*):", "\\1", table.features.raw)
  if(!identical(table.features, bc.features))
    stop("Features in the bicluster table differ from those annotated above. Files may corrupt")
  table.conditions.length <- unique(sapply(vector.split, length))-1L
  if(!identical(table.conditions.length, length(bc.conditions)))
    stop("Conditions in the bicluster table differ from those annotated above. Files may corrupt")

  #### expression table
  ## bic.exprs <- do.call(rbind, lapply(vector.split, function(x) as.integer(x[-1])))
  ## if(!identical(dim(bic.exprs), c(length(bc.features), length(bc.conditions))))
  ##  stop("The bicluster table is of wrong dimension. File may corrupt")

  bic.info <- list(Svalue=bc.Svalue,
                   features=bc.features,
                   conditions=bc.conditions)
  ##                ,exprs=bic.exprs)
  return(bic.info)
}

listCharToFactor <- function(list) {
  list.unique <- levels(factor(unlist(list)))
  list.index <- sapply(list, function(x) list.unique %in% x)
  return(list(factor=list.unique,
              indices=list.index))

}
isClassSanity <- function(input.list, class) {
  input.list.unique.class <- unique(sapply(input.list, class))
  identical(input.list.unique.class, class)
}


QUBICBicluster <- function(features, conditions, exprs, Svalue) {
  fdata <- data.frame(features=features)
  colnames(fdata)[1] <- RQUBIC_INPUT_FEATURENAME_SYMBOl
  fd <- new("AnnotatedDataFrame", fdata)
  pdata <- data.frame(conditions=conditions)
  colnames(pdata)[1] <- RQUBIC_INPUT_SAMPLENAME_SYMBOl 
  pd <- new("AnnotatedDataFrame", pdata)
  esObj <- new("QUBICBicluster",
               phenoData=pd, featureData=fd, exprs=exprs, Svalue=Svalue)
}

QUBICBiclusterSet <- function(qubicVersion, inputDataFile,
                              parameters, BCs, features, conditions) {
  parameterNames <- names(parameters)
  if(!setequal(parameterNames,g_QUBIC_PARAMETERS))
    stop("Missing or non-existing parameter(s) detected:\n",
         "Allowed parameter:", paste(g_QUBIC_PARAMETERS, collapse=","),"\n",
         "Input parameter:", paste(parameterNames, collapse=","))
  if(!identical(unique(sapply(BCs, class)), "QUBICBicluster"))
    stop("BCs must be a list of QUBICBicluster")
  res <- new("QUBICBiclusterSet",
             qubicVersion=qubicVersion,
             inputDataFile=inputDataFile,
             parameters=parameters,
             BCs=BCs,
             features=features,
             conditions=conditions)
  return(res)
}
  
parseQubicBlocks <- function(filename) {
  blocks <- scan(filename, what="character", sep="\n", quiet=TRUE)

  ## comment lines
  comment.lines <- grep("^#", blocks, value=TRUE)
  qubic.version <- gsub(".*version\\s*(([[:digit:]]|\\.)*)\\s*output", "\\1", comment.lines[1])
  qubic.datafile <- gsub(".*Datafile\\s*([[:print:]]*):.*", "\\1", comment.lines[2])
  qubic.datafile.type <- gsub(".*:\\s*(\\w*)\\stype$", "\\1", comment.lines[2])
  qubic.parameters <- gsub("#\\s*Parameters:\\s-k\\s(\\d*)\\s-f\\s([[:print:]]*)\\s-c\\s([[:print:]]*)\\s-o\\s(\\d*)\\s-q\\s([[:print:]]*)\\s-r\\s(\\d*)",
                           "\\1;\\2;\\3;\\4;\\5;\\6", comment.lines[3])
  qubic.parameters.split <- strsplit(qubic.parameters, ";")[[1]]
  qubic.para.k <- as.integer(qubic.parameters.split[1])
  qubic.para.f <- as.numeric(qubic.parameters.split[2])
  qubic.para.c <- as.numeric(qubic.parameters.split[3])
  qubic.para.o <- as.integer(qubic.parameters.split[4])
  qubic.para.q <- as.numeric(qubic.parameters.split[5])
  qubic.para.r <- as.integer(qubic.parameters.split[6])

  ## bc
  bc.heads.index <- grep("^BC", blocks, value=FALSE)
  bc.head.pattern <- "^BC\\d*\\sS=(\\d*)$"
  bc.Svalues <- as.integer(gsub(bc.head.pattern, "\\1", blocks[bc.heads.index]))
  bc.index <- seq(along=bc.heads.index)-1 ## note that the index of BCs starts with 0

  bc.nextheads.index <- c(bc.heads.index[-1]-1, length(blocks))
  bc.lines <- lapply(seq(along=bc.heads.index), function(i) seq(from=bc.heads.index[i],
                           to=bc.nextheads.index[i]))

  qubic.bics <- lapply(bc.lines, function(i) parseQubicBicluster(blocks[i]))
  qubic.length <- length(qubic.bics)
  qubics.bics.Svalues <- sapply(qubic.bics, function(x) x$Svalue)
  qubics.bics.features <- listCharToFactor(sapply(qubic.bics, function(x) x$features))
  qubics.bics.conditions <- listCharToFactor(sapply(qubic.bics, function(x) x$conditions))
  #### expression is not parsed, since can be access by subsetting
  ## qubics.bics.exprs <- lapply(qubic.bics, function(x) x$exprs)

  ## qubic.bicObjs <- lapply(seq(along=qubic.bics),
  ##                      function(x) QUBICBicluster(qubics.bics.features[[x]],
  ##                                                 qubics.bics.conditions[[x]],
  ##                                                 qubics.bics.exprs[[x]],
  ##                                                 qubics.bics.Svalues[[x]]))

  obj <- new("QUBICBiclusterSet",
             Parameters=list(
               Call=match.call(),
               Method="QUBIC",
               qubicVersion=qubic.version,
               inputDataFile=qubic.datafile,
               k=qubic.para.k,
               f=qubic.para.f,
               c=qubic.para.c,
               o=qubic.para.o,
               q=qubic.para.q,
               r=qubic.para.r),
             RowxNumber=qubics.bics.features[[2]],
             NumberxCol=t(qubics.bics.conditions[[2]]),
             Number=qubic.length,
             info=list(features=qubics.bics.features[[1]],
               conditions=qubics.bics.conditions[[1]],
               Svalues=qubics.bics.Svalues)
             )
  
  return(obj)
}

fileFound <- function(x) {
  if(!file.exists(x)) {
    warning(sprintf("File %s not found!\n", x))
    return(FALSE)
  }
  return(TRUE)
}

