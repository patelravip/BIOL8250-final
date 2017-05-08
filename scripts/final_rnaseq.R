## Start fresh!!
rm(list=ls())
cat("\014")

## Load packages
library("DESeq2")
library("tximport")

## Prepare vectors with sample information
analysisType = "rna_seq"
sampleNames = c("scer_rep1","scer_rep2","spar_rep1","spar_rep2")
sampleTypes = c("scer","scer","spar","spar")
sampleFiles = c(paste("data",sampleTypes,analysisType,sampleNames,"abundance.2.tsv", sep="/"))

## Read in files
txDat = tximport(sampleFiles, type="kallisto", txOut=TRUE)
coldata = data.frame(condition=sampleTypes)
rownames(coldata) = sampleNames

## DESeq2 object
dds = DESeqDataSetFromTximport(txDat, colData=coldata, design=~condition)
dex = DESeq(dds)
res = results(dex)
write.table(res, file="rnaseq_results.txt", sep="\t", quote=FALSE)

# ## Plot
# plotMA(res)
# # sizeFactorEst = estimateSizeFactors(dex)
# # dispEst = estimateDispersions(dex)
# plotDispEsts(dex, 1)