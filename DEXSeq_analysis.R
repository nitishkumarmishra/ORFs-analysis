# ###################################################################
# Author: Nitish Mishra
# Copyright (c) Nitish Mishra, 2022
# Email:  nitishimtech@gmail.com
# Date: 29 September, 2022
# Script Name: code_name.R
#####################################################################
# ##################  Script Description: ###########################
#
##################### Notes: ########################################
#
#####################################################################
################## SET WORKING DIRECTORY ############################
cat("SETTING WORKING DIRECTORY...\n\n", sep = "")
wd <- "C:/Users/nmishra/Dropbox/PC/Desktop/DEXSeq Analysis"
setwd(wd)
cat("WORKING DIRECTORY HAS BEEN SET TO: ", wd, sep = "")
#####################################################################
#setwd("C:/Users/nmishra/Dropbox/PC/Desktop/DEXSeq Analysis")
suppressMessages(suppressWarnings(library("DEXSeq")))
suppressMessages(suppressWarnings(library("BiocParallel")))

countFiles = list.files(getwd(), pattern="fb.txt$", full.names=TRUE)
countFiles = list.files(getwd(), pattern=".txt$", full.names=TRUE)
basename(countFiles)
flattenedFile = list.files(getwd(), pattern="gff$", full.names=TRUE)
flattenedFile


sampleTable = data.frame(
  row.names = gsub(".txt","", basename(countFiles)),
  condition = c("control", "control", "control",  
                "TGFb", "TGFb", "TGFb", "CX5461" , "CX5461", "CX5461"),
  libType = c( "paired-end", "paired-end", "paired-end", 
               "paired-end", "paired-end", "paired-end", "paired-end" , "paired-end", "paired-end") )



dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )

# Parallelization and large number of samples
# This part will work only in Linux environment
#BPPARAM = MultiCoreParam(4)

dxd = estimateSizeFactors( dxd )
#dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)


plotDispEsts( dxd )


# Testing for differential exon usage
dxd = testForDEU( dxd )
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
#dxd = testForDEU( dxd, BPPARAM=BPPARAM)
#dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
dxr1 = DEXSeqResults( dxd )
dxr1

mcols(dxr1)$description
table ( dxr1$padj < 0.1 )
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )
plotMA( dxr1, cex=0.8 )

sampleAnnotation(dxd)

## Save result in dataframe 
DEXSeq.results <- as.data.frame(dxr1)
DEXSeq.results <- DEXSeq.results[!is.na(DEXSeq.results$padj < 0.05), ]
DEXSeq.results.padj.0.05 <- DEXSeq.results[DEXSeq.results$padj < 0.05, ]
DEXSeq.results.padj.0.05 <- DEXSeq.results.padj.0.05[order(DEXSeq.results.padj.0.05$pvalue),]

plotDEXSeq( dxr1, "ENSMUSG00000000489", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr1, "ENSMUSG00000000489", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


plotDEXSeq( dxr1, "ENSG00000105379", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


save.image("DEXSeq_analysis.RData")
