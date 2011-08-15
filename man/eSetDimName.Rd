\name{eSetDimName}
\alias{eSetDimName}
\title{
  Get dimname from an eSet object
}
\description{
  This function is implemented to automatically validate
  and choose feature (sample) names from the user input. This function
  is exported for the purpose of easing other Bioconductor developers
  performing the similar job, and is not tended to be called by end-user
  directly. 
}
\usage{
eSetDimName(eset, input, type = c("feature", "sample"))
}
\arguments{
  \item{eset}{An object of \code{\linkS4class{eSet}} class, mostly an
    \code{\linkS4class{ExpressionSet}} class.}
  \item{input}{The user input, see details below}
  \item{type}{Either \sQuote{feature} or \sQuote{sample}, indicating
    which dimension should be determined}
}
\details{
  The input can be one of the following three possibilities:
  \itemize{
    \item Missing. Depending on the type, the results of calling
    \code{\link{featureNames}} (\dQuote{feature}) or
    \code{\link{sampleNames}} (\dQuote{sample}) on the \code{eset} object
    will be returned.
    \item A character string of length 1. Depending on the type, it is
    first to be machted to the column names of either
    \code{\link{fData}} or \code{\link{pData}} results of the
    \code{eset} object. If found, the values in that column are returned
    (coerced to characters if necessary). If not found, the function
    stops by raising an error.
    \item A character vector of the length equal to one of the two
    dimensions of the \emph{eset}. In this scenario, the function only
    validates the equality of the length, coerces the input into
    characters, and return them.
  }
  If none of the scenarios above was met, the function stops by raising
  an error.
}
\note{
  A special case arises if one of the dimensions of the \code{eset}
  object is \eqn{1}: In this case, the input value is interpreted as the
  new name and returned. No column name match will take place in this case.
}
\value{A vector of characters, the length of which determined by the
  dimension of the input object.}
\author{Jitao David Zhang <jitao_david.zhang@roche.com>}
\seealso{\code{\link{sampleNames}}, \code{\link{featureNames}},
  \code{\link{fData}}, \code{\link{pData}}

  \code{\link{writeQubicInputFile}} calls the function.
}
\examples{
data(sample.ExpressionSet)
sub.eset <- sample.ExpressionSet[1:3, 1:3]

## usage one:
eSetDimName(sub.eset, type="feature")
eSetDimName(sub.eset, type="sample")

## usage two
## "sex" is one column in the pData(sub.eset) 
eSetDimName(sub.eset, input="sex", type="sample")
\dontrun{eSetDimName(sub.eset, input="foo", type="sample")}

## usage three
eSetDimName(sub.eset, input=paste("Sample", 1:3), type="sample")
\dontrun{eSetDimName(sub.eset, input=paste("Sample", 1:4),
type="sample")}

## special case: dim equals to one
eSetDimName(sub.eset[,1], input="foo", type="sample")
}
