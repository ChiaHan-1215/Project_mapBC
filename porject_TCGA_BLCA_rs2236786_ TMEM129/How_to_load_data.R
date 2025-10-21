# Goal: download files from IMVigor210 package
# env: in Biowulf, R version 4.5.0 (2025-04-11)
# Data:10202025
# IMvigor210CoreBiologies INstall
# ref: http://research-pub.gene.com/IMvigor210CoreBiologies/#downloading-the-imvigor210corebiologies-package

library("biomaRt")
library("circlize")
library("ComplexHeatmap")
library("corrplot")
library("DESeq2")
library("dplyr")
library("DT")
library("edgeR")
library("ggplot2")
library("limma")
library("reshape2")
library("limma")
library("spatstat")
library("survival")
library("plyr")

#install.packages("DESeq_1.38.0.tar.gz", repos = NULL, type = "source")

library("DESeq")

#install.packages("IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL, type = "source")

library("IMvigor210CoreBiologies")


# Now we load data needed 
data(cds)

# counts 
cds@assayData[["counts"]]

# raw count
raw_genecount <- counts(cds)
# geneID list 
gene_info <- fData(cds)
# patient info
patient_info <- pData(cds)

# so in this folder: "/gpfs/gsfs12/users/leec20/R/rhel8/4.5/IMvigor210CoreBiologies/analysis"
# the Figure1.Rmd have code for normalized script


#https://www.sciencedirect.com/science/article/pii/S2059702923008463?via%3Dihub#abs0015
#Raw count data for the genes
#of interest were transformed to log2-normalized reads per
#million, and values for each gene were median centered
#across a representative reference clinical population

# load
data(cds)  
cds2 <- cds
data(fmone)  
fmi <- fmone



# normalize 
geneNames <- setNames(fData(cds2)$Symbol, 
                      as.character(rownames(fData(cds2))))

voomD <- filterNvoom(counts(cds2),
                     minSamples=ncol(counts(cds2))/10,
                     minCpm=0.25)

# the m in here is log2-counts per million (logCPM)
m <- voomD$E

# row-wise Z-score normalization per gene across samples.
m <- t(scale( t( m ),
              center=TRUE, 
              scale=TRUE)
)

# add signature scores to pData()
m2 <- m
rownames(m2) <- geneNames[rownames(m2)]



# https://support.bioconductor.org/p/40297/