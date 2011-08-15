library(rqubic)
source("test_utilities.R")

## parse QUBIC rules
rule.file <- getTestFile("sampleExpressionSet.rules")
rule.data <- parseQubicRules(rule.file)

## parse QUBIC chars
chars.file <- getTestFile("sampleExpressionSet.chars")
chars.data <- parseQubicChars(chars.file)

## parse QUBIC blocks
block.file <- getTestFile("sampleExpressionSet.blocks")
block.data <- parseQubicBlocks(block.file)

## QUBICBiclusterSet methods
length(block.data)
dim(block.data)
block.data.rxn <- RowxNumber(block.data)
block.data.nxc <- NumberxCol(block.data)
Svalue(block.data)
Svalue(block.data, 1:5)
Svalue(block.data, TRUE)

block.features <- features(block.data)
featureCount(block.data)
block.bcFeatures <- BCfeatures(block.data)
BCfeatures(block.data, 3:5)
BCfeatureCount(block.data)
BCfeatureCount(block.data, 3:5)

block.conds <- conditions(block.data)
conditionCount(block.data)
block.bcConditions <- BCconditions(block.data)
BCfeatures(block.data, 3:5)
BCconditionCount(block.data)
BCconditionCount(block.data, 3:5)

block.features <- featureNames(block.data)
block.conds <- sampleNames(block.data)

show(block.data)

## subset
block.data.sub <- block.data[1:3]

