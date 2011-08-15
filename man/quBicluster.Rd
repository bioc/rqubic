\name{quBicluster}
\alias{quBicluster}
\title{Qualitative Biclustering}
\description{
  The function takes seeds and quantileDiscretized ExpressionSet as input,
  biclusters the data and returns an object holding biclusters. Users
  may control the report number of clusters, tolerance of incoherent
  genes (or conditions), as well as the filtering of redundant clusters.
}
\usage{
quBicluster(seeds, eset, report.no = 100L, tolerance = 0.95, filter.proportion = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seeds}{An object of the S3-class \code{rqubicSeeds},
    representing seeds generated from the quantileDiscretized expression data}
  \item{eset}{Discretized expression data}
  \item{report.no}{Number of biclusters that should be
    reported. Detected biclusters are ranked by the S-score, which is
    defined by the product of gene counts and sample counts. They are
    ordered and the top ones are reported.}
  \item{tolerance}{Percentage of tolerated incoherent samples, 0.95 by default}
  \item{filter.proportion}{Proportion of a cluster, over which the
    cluster is considered as redudant. Each bicluster is compared to all
    better ranking biclusters, and the overlapping proportion is
    measured by the proportion of the product of overlapping samples and
    overlapping	genes, to the product samples and genes. If the
    proportion is larger than the given threshold, the block will be considered redundant and
  therefore not reported. Setting the threshold to 1 (default) does not
  perform any filtering.}
}
\details{
  The function calls a C routine to perform the biclustering. Currently
  the routine returns blocks with fewer samples specified by the minimum
  column number, due to the set of tolerance values. This might be changed in the fewer versions.
}
\value{
  An object of the \code{\link{QUBICBiclusterSet-class}}, holding all biclusters.
}
\references{
  Li et al. (2009) \emph{QUBIC: a qualitative biclustering algorithm for
  analyses of gene expression data} Nucleic Acids Research 37:e101
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\seealso{
\code{\link{quantileDiscretize}} and \code{\link{generateSeeds}}
}
\examples{
data(sample.ExpressionSet)
rqubic.example.discret <- quantileDiscretize(sample.ExpressionSet, rank=2L)

rqubic.example.sel.seeds <- generateSeeds(rqubic.example.discret, minColWidth=2L)

rqubic.example.blocks <- quBicluster(rqubic.example.sel.seeds,
    rqubic.example.discret,
    report.no=200L,
    filter.proportion=0.1)

## print features in each bicluster
BCfeatures(rqubic.example.blocks)
}
