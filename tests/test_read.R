library(rqubic)
source("test_utilities.R")

bf1 <- getTestFile("bic_1.out")
bfn <- getTestFile("bic_null.out")
bf3 <- getTestFile("bic_3.out")

bi1 <- readBiclusterResults(bf1)
bi3 <- readBiclusterResults(bf3)
class(bicn <- readBiclusterResults(bfn)) ## show method prints error

## combine
combineBiclusts(bi3, bi3)

## combine with empty bicluster objects
combineBiclusts(bi3, bi3, bicn)
