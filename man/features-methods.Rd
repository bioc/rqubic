\name{features-methods}
\docType{methods}

\alias{features}
\alias{featureCount}
\alias{conditions}
\alias{conditionCount}

\alias{BCfeatureCount}
\alias{BCfeatures}
\alias{BCconditionCount}
\alias{BCconditions}

\alias{features-methods}
\alias{conditions-methods}
\alias{featureCount-methods}
\alias{conditionCount-methods}

\alias{featureCount,QUBICBiclusterSet-method}
\alias{conditionCount,QUBICBiclusterSet-method}

\alias{featureNames,QUBICBiclusterSet-method}
\alias{sampleNames,QUBICBiclusterSet-method}

\alias{features,QUBICBiclusterSet-method}
\alias{featureCount,QUBICBicluster-method}
\alias{BCfeatures,QUBICBiclusterSet,missing-method}
\alias{BCfeatureCount,QUBICBiclusterSet,missing-method}
\alias{BCfeatures,QUBICBiclusterSet,ANY-method}
\alias{BCfeatureCount,QUBICBiclusterSet,ANY-method}

\alias{conditions,QUBICBiclusterSet-method}
\alias{conditionCount,QUBICBicluster-method}
\alias{BCconditions,QUBICBiclusterSet,missing-method}
\alias{BCconditionCount,QUBICBiclusterSet,missing-method}
\alias{BCconditions,QUBICBiclusterSet,ANY-method}
\alias{BCconditionCount,QUBICBiclusterSet,ANY-method}

\title{ Extract features and conditions }
\description{
  Generic function \code{features} and \code{conditions}, as well as auxillary
  count functions, are implemented for
  \code{\linkS4class{QUBICBiclusterSet}} objects.

  They can be used in one of the following forms:
  \enumerate{
    \item Used on a \code{\linkS4class{QUBICBiclusterSet}}, and without
    specifying index, \code{features} or \code{conditions} returns the unique and ordered
    features or conditions involved in at least one bicluster, and
    \code{featureCount} or \code{conditionCount} returns the length of
    repsective vectors. To get the feature/condition numbers in each
    bicluster of the set, use
    \code{\link{BCfeatureCount}}/\code{\link{BCconditionCount}}.
    \item Used on a \code{\linkS4class{QUBICBiclusterSet}} and provided
    index (indices), the features/conditions or their counts are
    returned for specified biclusters.
  }

  In addition, \code{\link{featureNames}} and \code{\link{sampleNames}}
  are of the same implementation as \code{\link{features}} and
  \code{\link{conditions}}.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "QUBICBicluster")}}{
  Information about all the biclusters.
}

\item{\code{signature(object = "QUBICBiclusterSet", index = "missing")}}{
  Information about all the biclusters in the set.
}

\item{\code{signature(object = "QUBICBiclusterSet", index = "ANY")}}{
  Information about selected biclusters in the set, the index can be
  integers or logical variables for subsetting.
}}
}
\references{
  Guojun Li, Qin Ma, Haibao Tang, Andrew H. Paternson and Ying Xu
  (2009) \emph{QUBIC: a qualitative biclustering algorithm for analyses
  of gene expression data}. Nucleic Acids Research, 37:e101
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\seealso{
  Methods are also documented in \code{\linkS4class{QUBICBiclusterSet}}.
}
\examples{
example.file <- system.file("extdata", "sampleExpressionSet.blocks", package="rqubic")

example.block <- parseQubicBlocks(example.file)

head(features(example.block))
featureCount(example.block)
head(conditions(example.block))
conditionCount(example.block)

BCfeatureCount(example.block)
BCfeatures(example.block)[1:2]
BCconditionCount(example.block)
BCconditions(example.block)[1:2]

head(featureNames(example.block))
head(sampleNames(example.block))
}
\keyword{methods}
