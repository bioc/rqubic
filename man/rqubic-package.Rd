\name{rqubic-package}
\alias{rqubic-package}
\alias{rqubic}
\docType{package}
\title{
  Qualitative biclustering algorithm for expression data analysis
}
\description{
  QUBIC is a qualitative biclustering algorithm for high-throughput
  expression data analysis. \code{rqubic} package implements this
  algorithm in R, partly with the codes contributed by Haibao Tang and Qin Ma (version 0.23 released without
  any limitation).

  The \code{rqubic} package also provides parsers for the command
  line tool of qubic written in C.
}
\details{
\tabular{ll}{
Package: \tab rqubic\cr
Type: \tab Package\cr
Version: \tab 1.5\cr
Date: \tab 2011-04-11\cr
License: \tab LGPL-2 \cr
LazyLoad: \tab yes\cr
}

Part of the source code in C is modified from the source code of the QUBIC
command line tool (in C) provided by Haibao Tang and Qin Ma
\email{maqin@csbl.bmb.uga.edu}, downloaded from
\url{http://csbl.bmb.uga.edu/~maqin/bicluster/} on 01.03.2011, version 0.23.

Source code of QUBIC also uses open-source data structure library
codes. See the README file included in the QUBIC command line tool source.
}
\author{
Jitao David Zhang <jitao_david.zhang@roche.com>, Laura Badi and Martin Ebeling
Maintainer: Jitao David Zhang <jitao_david.zhang@roche.com>
}
\references{
  Li et al. (2009) \emph{QUBIC: a qualitative biclustering algorithm for
  analyses of gene expression data} Nucleic Acids Research 37:e101
}
\keyword{ package }
%%\examples{
%%
%%}
