#!/usr/bin/Rscript

# Script information ------------------------------------------------------

# title: Illumina HT-12 array gene expression data processing 
# author: Walter Muskovic
# date: 2020-08-11
# description: Process gene expression data for 35 DIPG and 10 normal brain samples generated using Illumina's Human HT-12 WG-DASL V4.0 expression beadchip technology



# Import libraries --------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lumi)
  library(lumiHumanIDMapping)
})



# Download data -----------------------------------------------------------

download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50021/suppl/GSE50021_non-normalized.txt.gz",
              destfile = "data/GSE50021_non-normalized.txt.gz")



# Process data ------------------------------------------------------------

# Read Bead Studio output file and create a LumiBatch object
expression_data <- lumiR.batch('data/GSE50021_non-normalized.txt.gz', lib.mapping = "lumiHumanIDMapping")

# summary of data
expression_data

# Perform quantile normalization (No variance stabilizing or log transform)
expression_data_norm <- lumiExpresso(expression_data, variance.stabilize=FALSE, normalize.param = list(method='quantile'))

# Plot expression value distributions before and after normalisation
png(filename = "figures/normalisation.png", width = 1000, height=800)
par(mfrow=c(2,1))
boxplot(exprs(expression_data), main="Before normalisation", las=2)
boxplot(exprs(expression_data_norm), main="After normalisation")
par(mfrow=c(1,1))
dev.off()



# Write out data ----------------------------------------------------------

# Get the expression data
expression_matrix <- exprs(expression_data_norm)

# Add gene names
nuIDs <- featureNames(expression_data_norm)
mappingInfo <- getNuIDMappingInfo(nuIDs, lib.mapping='lumiHumanIDMapping')
gene_name <- mappingInfo[,"Symbol"][match(row.names(expression_matrix), row.names(mappingInfo))] %>% as.character
expression_matrix <- cbind(data.frame(gene_name), round(data.frame(expression_matrix), 2))

# Write out as a tab-delimited file
write.table(expression_matrix, "data/expression_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
