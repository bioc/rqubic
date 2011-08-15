##======================================================================##
## Data Structures of rqubic package
##
## Author: Jitao David Zhang <jitao_david.zhang@roche.com>, March 15, 2011
##======================================================================##
## global variables
g_QUBIC_PARAMETERS <-  c("k", "f", "c", "o", "q", "r")
RQUBIC_INPUT_FEATURENAME_SYMBOl <- "RQUBIC_INPUT_FEAT"
RQUBIC_INPUT_SAMPLENAME_SYMBOl <- "RQUBIC_INPUT_SAMPLE"

setClass("QUBICBiclusterSet", contains="Biclust")
