\name{parseQubicBlocks}
\alias{parseQubicRules}
\alias{parseQubicChars}
\alias{parseQubicBlocks}
\title{
  Parse QUBIC Command Line Tool Output Files
}
\description{
  These functions parse output files of the QUBIC command line tool
  developed by Ma et al. 
}
\usage{
parseQubicRules(filename)
parseQubicChars(file, check.names=FALSE, ...)
parseQubicBlocks(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{Input filename}
  \item{file}{Input filename}
  \item{check.names}{logical, should the column names be checked?}
  \item{...}{other parameters passed to the \code{\link{read.csv}} function}
}
\details{
  Parse QUBIC Command Line Tool Output Files
}
\value{
  \code{parseQubicRules} and \code{parseQubicChars} both return a data
  frame.

  \code{parseQubicBlocks} returns an instance of
  \code{\linkS4class{QUBICBiclusterSet}} class.
}
\references{
  \url{http://csbl.bmb.uga.edu/~maqin/bicluster/}
}
\author{
  Jitao David Zhang <jitao_david.zhang@roche.com>
}
\examples{
getRqubicFile <- function(filename) system.file("extdata", filename, package="rqubic")

## parse QUBIC rules
rule.file <- getRqubicFile("sampleExpressionSet.rules")
rqubic.sample.rule <- parseQubicRules(rule.file)

## parse QUBIC chars
chars.file <- getRqubicFile("sampleExpressionSet.chars")
rqubic.sample.chars <- parseQubicChars(chars.file)

## parse QUBIC blocks
block.file <- getRqubicFile("sampleExpressionSet.blocks")
rqubic.sample.data <- parseQubicBlocks(block.file)
}
