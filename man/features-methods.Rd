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

\alias{featureCount,Biclust-method}
\alias{conditionCount,Biclust-method}

\alias{featureNames,Biclust-method}
\alias{sampleNames,Biclust-method}

\alias{features,Biclust-method}
%%\alias{featureCount,QUBICBicluster-method}
\alias{BCfeatures,Biclust,missing-method}
\alias{BCfeatureCount,Biclust,missing-method}
\alias{BCfeatures,Biclust,ANY-method}
\alias{BCfeatureCount,Biclust,ANY-method}

\alias{conditions,Biclust-method}
%%\alias{conditionCount,QUBICBicluster-method}
\alias{BCconditions,Biclust,missing-method}
\alias{BCconditionCount,Biclust,missing-method}
\alias{BCconditions,Biclust,ANY-method}
\alias{BCconditionCount,Biclust,ANY-method}

\title{ Extract features and conditions }
\description{
  Generic function \code{features} and \code{conditions}, as well as auxillary
  count functions, are implemented for
  \code{\linkS4class{Biclust}} objects.

  They can be used in one of the following forms:
  \enumerate{
    \item Used on a \code{\linkS4class{Biclust}}, and without
    specifying index, \code{features} or \code{conditions} returns the unique and ordered
    features or conditions involved in at least one bicluster, and
    \code{featureCount} or \code{conditionCount} returns the length of
    repsective vectors. To get the feature/condition numbers in each
    bicluster of the set, use
    \code{\link{BCfeatureCount}}/\code{\link{BCconditionCount}}.
    \item Used on a \code{\linkS4class{Biclust}} and provided
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

\item{\code{signature(object = "Biclust", index = "missing")}}{
  Information about all the biclusters in the set.
}

\item{\code{signature(object = "Biclust", index = "ANY")}}{
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
  See other methods implemented for the \code{\linkS4class{Biclust}}
  class in the \code{biclust} package. And the methods specific for
  \code{\linkS4class{QUBICBiclusterSet}}.
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
