\name{fcFilter}
\alias{fcFilter}
\alias{fcFilter,Biclust-method}
\title{
  Feature-Condition Filter
}
\description{
  Filter Biclusters by feature and concition counts. Biclusters with
  fewer features or conditions than specified thresholds are removed.
}
\usage{
fcFilter(object, ...)
}
\arguments{
  \item{object}{A \code{Biclust} object}
  \item{\dots}{Two parameters are accepted: \code{feat.min} and
    \code{cond.min}. They indicate the minimum number of features
    (conditions) in biclusters, which are to be accepted. It is allowed to only specify one threshold. If
    both are specified, only biclusters fulfilling both criterion are accepted.}

}
\value{
  A \code{Biclust} object.
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
data(sample.ExpressionSet)
rqubic.example.discret <- quantileDiscretize(sample.ExpressionSet, rank=2L)

rqubic.example.sel.seeds <- generateSeeds(rqubic.example.discret, minColWidth=2L)

rqubic.example.blocks <- quBicluster(rqubic.example.sel.seeds,
    rqubic.example.discret,
    report.no=200L,
    filter.proportion=0.1)

print(rqubic.example.blocks)
print(fcFilter(rqubic.example.blocks,feat.min=10))
print(fcFilter(rqubic.example.blocks,cond.min=2))
print(fcFilter(rqubic.example.blocks,feat.min=10, cond.min=2))
}
