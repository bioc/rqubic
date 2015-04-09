\name{readBiclusterResults}
\alias{readBiclusterResults}
\title{Import bicluster results from plain text files}
\description{
  This function complements the functionality of
  \code{\link{writeBiclusterResults}} in the \code{biclust} package. It
  constructs a \code{\linkS4class{Biclust}} object from a plain text file.
}
\usage{
readBiclusterResults(filename, featureNames, sampleNames, delimiter = ";", ...)
}
\arguments{
  \item{filename}{Character, name of the file storing biclustering information}
  \item{featureNames}{Optional character vector, feature names of the
    underlying expression dataset. See details.}
  \item{sampleNames}{Optional character vector, sample names of the
    underlying expression dataset. See details.}
  \item{delimiter}{Character used to separate features, samples and
    counts of them.}
  \item{\dots}{Other parameters passed to the \code{\link{readLines}} function.}
}
\details{
  Currently output files written by the \code{writeBiclusterResults}
  function does not contain original feature names or sample names in
  the expression dataset from which biclusters were mined. The
  \code{featureNames} and \code{sampleNames} allow to use this
  information to construct a \code{\linkS4class{Biclust}} object that
  has the same dimension as the original expression dataset.
}
\value{
  A \code{\linkS4class{Biclust}} object 
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\seealso{
  In case several biclustering algorithms were applied to the same
  expression dataset, they can be combined with
  \code{\link{combineBiclusts}} once the results were read from plain texts.
}
\examples{
temp <- tempfile()
library(biclust)
data(BicatYeast, package="biclust")
res <- biclust(BicatYeast, method=BCBimax(), number=5)
writeBiclusterResults(temp, res,"CC with delta 1.5",dimnames(BicatYeast)[[1]],dimnames(BicatYeast)[[2]], delimiter=";")
res.back <- readBiclusterResults(temp, delimiter=";")
res.full.back <- readBiclusterResults(temp,featureNames=rownames(BicatYeast), sampleNames=colnames(BicatYeast),delimiter=";")
}

