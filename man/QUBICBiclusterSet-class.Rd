\name{QUBICBiclusterSet-class}
\Rdversion{1.1}
\docType{class}
\alias{QUBICBiclusterSet-class}
\alias{Svalue,QUBICBiclusterSet,missing-method}
\alias{Svalue,QUBICBiclusterSet,ANY-method}
\alias{Svalue,eSet,missing-method}
\alias{Svalue,matrix,missing-method}
\alias{[,Biclust,ANY,missing,missing-method}

\alias{BCcount,Biclust-method}
\alias{info,Biclust,missing-method}
\alias{info,Biclust,ANY-method}
\alias{parameter,Biclust,missing-method}
\alias{parameter,Biclust,character-method}
\alias{RowxNumber,Biclust-method}
\alias{NumberxCol,Biclust-method}
\alias{show,QUBICBiclusterSet-method}

%% alias generics
\alias{BCcount}
\alias{Svalue}
\alias{info}
\alias{parameter}
\alias{RowxNumber}
\alias{NumberxCol}

\title{Class "QUBICBiclusterSet"}
\description{
  Object representing a set of biclusters identified by the QUBIC
  algorithm. The class structure inherits the
  \code{\linkS4class{Biclust}} class in the \code{biclust} package.
}
\section{Objects from the Class}{
  Created by functions parsing the output files of QUBIC command line
  tool, or functions calling QUBIC algorithm implementations in R.

  Not intended to be created manually by end-users. However, interested
  users are invited to review the source code or use the
  \code{\link{showClass}} method to view the construction of the class.
}
\section{Slots}{
  \describe{
    See the class structure of \code{\linkS4class{Biclust}}. The slots
    \code{Parameter} and \code{Info} have been filled with lists
    releveant to the QUBIC algorithm, and all items should be accessed
    by S4-methods to make sure the consistency.
  }
}
\section{Methods}{
  \describe{
    \item{Svalue}{\code{signature(object = "QUBICBiclusterSet", index =
	"missing")}: Return S values of QUBIC biclusters as a vector}
    \item{Svalue}{\code{signature(object = "QUBICBiclusterSet", index = "numeric")}: S
      values of specified bicluster(s) are returned. Index is one or a
      vector of integers. Non-integers will be coereced.}
    \item{[}{\code{signature(x = "QUBICBiclusterSet", i = "ANY", j =
	"missing", drop = "missing")}: Returning a
      subset of the current QUBICBiclusterSet.}
    \item{parameter}{\code{signature(object = "Biclust", index
	= "character")}: return an input parameter specified by the
      parameter name}
    \item{parameter}{\code{signature(object = "Biclust", index
	= "missing")}: return a list of input parameters used by the
      biclustering algorithm, for example QUBIC}
    \item{info}{\code{signature(object = "Biclust", index
	= "ANY")}: return information of the biclusters. For end-users,
      specific information accessors should be preferred, for example
      \code{\link{features}}, \code{\link{conditions}} and \code{Svalue}}
    \item{info}{\code{signature(object = "Biclust", index
	= "missing")}: return all information of the biclusters in a list. For end-users,
      specific information accessors should be preferred, for example
      \code{\link{features}}, \code{\link{conditions}} and \code{Svalue}}
    \item{show}{\code{signature(object = "QUBICBiclusterSet")}: showing method }
  }
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
  \code{\linkS4class{Biclust}} is the basic block accomodating
  biclusters identified by the QUBIC algorithm.
}
\examples{
showClass("QUBICBiclusterSet")
}
\keyword{classes}
