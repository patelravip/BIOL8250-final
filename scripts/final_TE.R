## Start fresh!!
rm(list=ls())
cat("\014")

## Load packages
library("DESeq2")
library("tximport")

## Prepare vectors with sample information
sampleNames = c("scer_rep1","scer_rep2","spar_rep1","spar_rep2")
rnaSampleNames = paste(sampleNames,"_rna", sep="")
riboSampleNames = paste(sampleNames,"_ribo", sep="")
sampleTypes = c("scer","scer","spar","spar")
rnaSeqSampleFiles = c(paste("../data",sampleTypes,"rna_seq",sampleNames,"abundance.2.tsv", sep="/"))
riboSeqSampleFiles = c(paste("../data",sampleTypes,"ribo_seq",sampleNames,"abundance.2.tsv", sep="/"))

## Concatenate appropriate lists
sampleFiles = c(rnaSeqSampleFiles, riboSeqSampleFiles)
assayTypes = c(rep("rna",4), rep("ribo",4))

## Read in files
txDat = tximport(sampleFiles, type="kallisto", txOut=TRUE)
coldata = data.frame(condition=rep(sampleTypes,2), assay=assayTypes)
rownames(coldata) = c(rnaSampleNames,riboSampleNames)

## DESeq2 object for differential translational efficiency test
dds2 = DESeqDataSetFromTximport(txDat, colData=coldata, design=~assay+condition+assay:condition)
dds2$assay = relevel(dds2$assay, ref="rna")
dex2 = DESeq(dds2, test="LRT", reduced=~assay+condition)
res2 = results(dex2)
write.table(res2, file="../tables/TXT/te_table.txt", sep="\t", quote=FALSE)

# ## Plot
# plotMA(res, ylim=c(-3,3), xlim=c(1e-1,1e5))
# # sizeFactorEst = estimateSizeFactors(dex)
# # dispEst = estimateDispersions(dex)
# plotDispEsts(dex, 1)