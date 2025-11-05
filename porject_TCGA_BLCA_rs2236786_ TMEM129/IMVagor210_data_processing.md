### Goal: Expore the IMvagor210 package
#### env: in Biowulf, R version 4.5.0 (2025-04-11)
#### Data:10202025


#### IMvigor210CoreBiologies Install step:

ref: 

http://research-pub.gene.com/IMvigor210CoreBiologies/#downloading-the-imvigor210corebiologies-package

https://github.com/SiYangming/IMvigor210CoreBiologies

**Step1: create Rstduio env on Biowulf sintereactive mode**
As the biowulf have build-in packages required for IMvagor. 

We need to download the raw package file of DESeq (https://github.com/SiYangming/DESeq/releases/tag/v1.39.0) and IMvagor210 (http://research-pub.gene.com/IMvigor210CoreBiologies/packageVersions/) to install locally.

```R
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

# insatll DEseq 
#install.packages("./fodler_to_tools/DESeq_1.38.0.tar.gz", repos = NULL, type = "source")

library("DESeq")

# Install main package 
#install.packages("./fodler_to_tools/IMvigor210CoreBiologies_1.0.0.tar.gz", repos = NULL, type = "source")
library("IMvigor210CoreBiologies")

```

Once all library loaded successfully, now we can load data 

```R
# Now we load data needed 
data(cds)

# raw gene count
raw_genecount <- counts(cds)
# geneID list 
gene_info <- fData(cds)
# patient info
patient_info <- pData(cds)

```

in this folder: `/gpfs/gsfs12/users/leec20/R/rhel8/4.5/IMvigor210CoreBiologies/analysis` contains scripts for generating figures for thier paper.

In `Figure1.Rmd` have code for normalized their raw count to z-score value using Voom package.

```R

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

```

To match our analysis methods, thanks Oscar proveied the other normalized way and to transform to z-score.

the first is to get normalized gene count using Deseq2

```R
# using Deseq2 to output normlaized count
library(DESeq2)
cts <- as.matrix(raw_genecount)
# Optional: ensure integers
storage.mode(cts) <- "integer"
coldata <- DataFrame(row.names = colnames(cts))  # zero-column DataFrame with sample IDs

colData <- data.frame(
  sample_id = colnames(cts),
  row.names = colnames(cts)
)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ 1) # ~ 1 means no variables in the design


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)                 # or DESeq(dds) if you'll do DE later
norm_counts <- counts(dds, normalized = TRUE) 
norm_counts
rownames(norm_counts) <- geneNames[rownames(norm_counts)]
norm_counts <- as.data.frame(norm_counts)

norm_counts$gene_symbol <- rownames(norm_counts)
norm_counts <- norm_counts[,c(349,1:348)]


#####################Do z score on Oscar's code ##############################

norm_counts <- norm_counts[,-1]
str(norm_counts)

# The file format should look like this:

#TCGA.2F.A9KO TCGA.2F.A9KP TCGA.2F.A9KQ TCGA.2F.A9KR TCGA.2F.A9KT
# TACC3         4.9714       5.7981       4.9421       7.0476       6.9486
# FGFR3         6.0588       7.7760       7.4456       8.8113       8.9363
# SLBP          5.3292       6.1125       5.6256       5.8066       6.8386
# TMEM129       3.7971       5.4871       5.4430       6.0605       5.9451
# FAM53A        0.2522       1.7702       0.2400       0.8246       1.8036

# Need to matched with data format 

head(norm_counts[1:5,1:5])

qnorm.DeseqNC.set <- t(apply(norm_counts, 1, rank, ties.method = "average"))
qnorm.DeseqNC.set <- qnorm(qnorm.DeseqNC.set / (ncol(qnorm.DeseqNC.set) + 1))
qnorm.DeseqNC.set[1:5,1:5]
qnorm.DeseqNC.set <- as.data.frame(qnorm.DeseqNC.set)
qnorm.DeseqNC.set$Gene_symbol <- rownames(qnorm.DeseqNC.set)
qnorm.DeseqNC.set <- qnorm.DeseqNC.set[,c(349,1:348)]
#write.table(qnorm.DeseqNC.set,'IMvigor_gene_zscore_from_Deseq2_normalizedcount_usingOscarMethod.csv',col.names = T,row.names = F,quote = F,sep = ',')

```

Now we get the z-score of Deseq2 normalized gene count! 
