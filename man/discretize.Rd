\name{quantileDiscretize}
\alias{quantileDiscretize}
\alias{quantileDiscretize-methods}
\alias{quantileDiscretize,matrix-method}
\alias{quantileDiscretize,eSet-method}
\title{
  Discretize expression matrix for qualitative biclustering
}
\description{
  Performs recursive quantilizations on gene expression data across
  samples, to quantileDiscretize gene expression matrix. The quantile parameter
  \code{q} determines the estimated proportion of differentially
  expressed genes (\emph{2q} as for both up- and down-regulatons). The
  rank parameter \code{r} determines how many discrete levels should
  differentially expressed genes (or outliers) have. See details below.
}
\usage{
quantileDiscretize(x, ...)
}
\arguments{
  \item{x}{It can be an object of the \code{\linkS4class{eSet}} class or
    inheriting it. The most commonly used form is an
    \code{linkS4class{ExpressionSet}} class. Alternatively, it can be a
    numeric matrix.}
  \item{...}{Currently, the \dots accepts two parameter: \code{q} and
    \code{rank}, explained below.}
  \itemize{
    \item{q}{Estimated proportion of conditions where gene is up- or
      down-regulated, value between \eqn{(0,0.5)}, default value is set to 0.06. By specifying \code{q} one
      estimates that in \code{2q} of all conditions, the expression value
      of a gene is considered as outlier.}
    \item{rank}{Ranks (levels) of outliers, a positive integer, default is 1L. By default,
      all conditions get one label for each gene in \eqn{{-1, 0, 1}},
      representing down expression, not changing and high expression
      respectively. In case \eqn{rank>1}, the outliers are further divided into
      \emph{rank} levels by applying recursive quantilization with equal
      intervals.}
  }
}
\details{
  Parameter \code{q} corresponds to the command line option \code{-q}
  in the QUBIC command line tool, and the \code{rank} option corresponds to
  \code{-r}.

  For each gene, the algorithm applies quantile discretization first to
  divide conditions into negative (lower), un-changed and positive (higher) expressions. Negative
  and positive expressed conditions are considered as \emph{outliers}. For
  outliers in each direction, the algorithm tries to further quantileDiscretize
  the expression values in case \eqn{rank>1}.

  This second discretization step is performed by dividing the sorted
  outliers into \eqn{rank} tandom groups with equal conditions. A label
  is assigned to each of these tandom groups, in the following order:
  \deqn{-1, -2, \ldots, -rank}
  for outliers with negative expression, from the \emph{most negative group} to
  the \emph{least negative group} (not the other way around!).

  Similarly, for positive outliers, labels in the order of
  \deqn{rank, rank-1, \ldots, 1}
  are assigned to tandom groups from the \emph{least positive group} to
  the \emph{most positive group}.

  That is, signs of labels indicate the direction of gene expression
  change, and the absolute value represents the quantileDiscretized \emph{rank}
  in the outliers.
}
\value{
  An object of the same class as the input parameter, with the
  \code{exprs} slot replaced by the quantileDiscretized matrix, which is a
  matrix of integer.
}
\references{
  Li et al. (2009) \emph{QUBIC: a qualitative biclustering algorithm for
  analyses of gene expression data} Nucleic Acids Research 37:e101
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\note{
  Note that the resulting discrete matrix of this implementation can be
  slighly different from the one used by the QUBIC command line
  tool.

  The main reason for this is the internal data type: while QUBIC
  uses \code{float} to represent expression matrix, we use \code{double}
  to represent the matrix.

  It has the advantages of interfacing to R,
  having higher precision and avoiding errors caused by floating
  presentation. It is implemented with potential larger costs of memory,
  however for test data sets (for example the ALL dataset with more than
  120 samples and 12000 genes) the peak memory use (<100M) as well as
  the execution time (CPU time 0.028s) are well under control.

  The differentially is especially often observed when there are many
  tied values. These cases however are very rare cases and we assume
  they will not affect the results to a large extent.
}

\seealso{
  \code{\link{parseQubicChars}} parses the quantileDiscretized matrix by the
  QUBIC command line tool into a data frame.
}
\examples{
library(Biobase)
data(sample.ExpressionSet, package="Biobase")
sample.disc <- quantileDiscretize(sample.ExpressionSet)
exprs(sample.disc)[1:6, 1:6]

## Equivalent to pass a numeric matrix
sample.mat.disc <- quantileDiscretize(exprs(sample.ExpressionSet))
sample.mat.disc[1:6, 1:6]
\dontrun{identical(exprs(sample.disc),sample.mat.disc)}

## with multiple ranks
sample.rank3 <- quantileDiscretize(sample.ExpressionSet, rank=3)
exprs(sample.rank3)[1:6, 1:6]
}
