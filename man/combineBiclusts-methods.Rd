\name{combineBiclusts-methods}
\docType{methods}
\alias{combineBiclusts}
\alias{combineBiclusts-methods}
\alias{combineBiclusts,Biclust,Biclust-method}
\alias{combineBiclusts,QUBICBiclusterSet,QUBICBiclusterSet-method}
\alias{combineBiclusts,list,missing-method}
\title{ Combine two or more Biclust objects into one}
\description{
  Combining two or more \code{Biclust} objects into one object. These
  objects must be of same dimension, namely have same numbers of
  features and samples, although the numbers of biclusters do not have
  to be the same (and usually are not).
}
\section{Methods}{
\describe{

\item{\code{signature(x = "Biclust", y = "Biclust")}}{Method for any
  Biclust objects}

\item{\code{signature(x = "QUBICBiclusterSet", y =
    "QUBICBiclusterSet")}}{Method for \code{QUBICBiclusterSet}
  only. Besides combining biclusters, they will also combine parameters
  and information stored in each \code{QUBICBiclusterSet} into the
  returning object.}

\item{\code{signature(x = "list", y =
    "missing")}}{Method for a list of \code{Biclust} or
  \code{QUBICBiclusterSet} objects.}

}}
\keyword{methods}
\examples{
data(sample.ExpressionSet)
re1.discret <- quantileDiscretize(sample.ExpressionSet, rank=1L)
re1.sel.seeds <- generateSeeds(re1.discret, minColWidth=2L)
re1.blocks <- quBicluster(re1.sel.seeds,
    re1.discret,
    report.no=50L,
    filter.proportion=0.1)

re2.discret <- quantileDiscretize(sample.ExpressionSet, rank=2L)
re2.sel.seeds <- generateSeeds(re2.discret, minColWidth=2L)
re2.blocks <- quBicluster(re2.sel.seeds,
    re2.discret,
    report.no=50L,
    filter.proportion=0.1)

re3.discret <- quantileDiscretize(sample.ExpressionSet, rank=3L)
re3.sel.seeds <- generateSeeds(re3.discret, minColWidth=2L)
re3.blocks <- quBicluster(re2.sel.seeds,
    re2.discret,
    report.no=50L,
    filter.proportion=0.1)

re12.blocks <- combineBiclusts(re1.blocks, re2.blocks)
re123.blocks <- combineBiclusts(re1.blocks, re2.blocks, re3.blocks)
re123.list.blocks <- combineBiclusts(list(re1.blocks, re2.blocks,
re3.blocks))
stopifnot(identical(re123.blocks, re123.list.blocks))
}
