library(rqubic)
source("test_utilities.R")

## quantileDiscretize
data(sample.ExpressionSet)
exp.eset <- sample.ExpressionSet
fData(exp.eset) <- data.frame(label=rownames(exprs(exp.eset)),
                              row.names=rownames(exprs(exp.eset)))

## test for C code
## note that to let QUBIC use all the samples, the first cell of the rectangular table must be a "o"!
dis.rank <- 1L
exp.discret <- quantileDiscretize(exp.eset, rank=dis.rank)
stopifnot(identical(exprs(exp.discret), quantileDiscretize(exprs(exp.eset), rank=dis.rank)))

## Note that the return value of quantileDiscretize is not exactly the same as the version implemented in the C tool of QUBIC
## The reason is that in this implementation we use data type "double" to represent expression value, which has higher precision
##   than the "float" type used in the QUBIC tool

## The following codes visualizes the differences. They mostly appear where there are many tied values.
if(FALSE) {
  chars.mat <- as.matrix(chars.data[,-1])
  diff.max <- chars.mat - exprs(exp.discret)
  diff.max.rowMan <- apply(diff.max, 1, function(x) sum(x!=0))
  cbind(qubic=chars.mat[134,], rqubic=exprs(exp.discret)[134,])
  if(require(RColorBrewer))
    image(t(diff.max), col=brewer.pal(5, "RdBu"))
}

## Generate seeds
exp.sel <- 1:nrow(exp.discret)
exp.sel.mincol <- 2L
exp.discret.sel <- exp.discret[exp.sel,]
exp.sel.seeds <- generateSeeds(exp.discret.sel, minColWidth=exp.sel.mincol)
exp.sel.seeds

## writeQubicInputFile(exp.discret.sel, file="testCase.csv", sampleNames=1:ncol(exp.discret.sel))

## TODO
## quBicluster: a list of
exp.blocks <- quBicluster(exp.sel.seeds,exp.discret.sel,
                          report.no=100L,
                          filter.proportion=0.1)
exp.blocks
exp.blocks[1]
parameter(exp.blocks)
features(exp.blocks)
BCfeatures(exp.blocks)
conditions(exp.blocks)
BCconditionCount(exp.blocks)
