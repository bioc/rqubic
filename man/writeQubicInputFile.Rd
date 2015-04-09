\name{writeQubicInputFile}
\alias{writeQubicInputFile}
\title{
  Write an ExpressionSet object into the file format required by the
  QUBIC command line tool
}
\description{
  The QUBIC commmand line tool (developed by Ma et al.) requires a
  tab-limited data matrix as input file, with some special requirements
  (see details below). This function takes an object of ExpressionSet
  and outputs the file.
}
\usage{
writeQubicInputFile(x, file = "", featureNames, sampleNames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object inheriting the \code{\linkS4class{eSet}} class,
    most commonly an \code{\linkS4class{ExpressionSet}} object,
    representing expression data of features across samples.}
  \item{file}{Filename to output, or a connection to write to (e.g. \code{stdout()}).}
  \item{featureNames}{Specifies the feature names. It can be left blank,
  in which case the result of calling \code{\link{featureNames}}
  on \code{x} will be used. Alternatively, it can be one character
  string, specifying which column in the \code{\link{fData}} should be
  used. The third possibility, it can be a vector of characters, with
  the same length as features in the object. In the last option, all
  other types will be converted to characters.}
  \item{sampleNames}{Specifies the sample names. It can be left blank,
  in which case the result of calling \code{\link{sampleNames}}
  on \code{x} will be used. Alternatively, it can be one character
  string, specifying which column in the \code{\link{pData}} should be
  used. The third possibility, it can be a vector of characters, with
  the same length as features in the object. In the last option, all
  other types will be converted to characters.}
}
\details{
  The description of the data
  format can be checked by running the QUBIC tool in the command line
  mode, with the option \emph{-h} (for \emph{help}). A special
  requirement, which makes it different from the results of the
  \code{\link{write.table}} function in R, is that before the sample
  names (column names), an \dQuote{o} must be added.
}
\value{
  No visible value will be returned, the function is called for its side effect.
}
\author{Jitao David Zhang <jitao_david.zhang@roche.com>}
\seealso{\code{\link{eSetDimName}}, \code{\link{write.table}}}
\examples{
tmpfile <- tempfile()
data(sample.ExpressionSet, package="Biobase")
sub.eset <- sample.ExpressionSet[1:3, 1:3]

## write to standard output
writeQubicInputFile(sub.eset)

## write to a temporary file
writeQubicInputFile(sub.eset, tmpfile)
head(readLines(tmpfile))

## specify names with one column name in fData/pData
writeQubicInputFile(sub.eset, file="", sampleNames="sex")

## alternatively specifiy names manually
writeQubicInputFile(sub.eset, file="",sampleNames=paste("Sample", 1:3))
}

