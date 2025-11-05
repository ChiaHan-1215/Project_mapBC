## Goal: Extract the TMEM129 exon expression value from GSE87304 using GEOquary and other tools

**The script is still under develope and need to recheck some value 

For one of the paper method they mentioned: 
```
we obtained the data from GSE87304 using the GEOquery R package (v.2.54.1), reprocessed the CEL files using the frma R package (v1.40.0),  annotated using R package ensembledb (c. 2.12.1, Ensemble 99) and applied a log2 transformation.
```

We need to download the raw CEL file and reprocess it to get the exon value.

Download `GSE87304_RAW.tar` from their GSE link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87304

Unzip and extracted all the xxx.cel.gz file 

some of the cel file can not reanalysis to get exon level expression, will try to fogure out ohter way?

```R
library(GEOquery)
library(DECIPHER)
library(pdInfoBuilder)
library(oligo)
library(affxparser)
library(frma)
library(huex.1.0.st.v2frmavecs)

# set folder to cel files 
setwd('project_TMEM129_TCGA_BLCA/CEL_stuff/')

chip_types <- sapply(list.celfiles(".", full.names=TRUE), function(f)
       tryCatch(readCelHeader(f)$chiptype, error=function(e) NA))

# so there are 7 data not fit for exon-level
table(chip_types)

# load only fitted data 
huex_files <- names(chip_types)[chip_types == "HuEx-1_0-st-v2"]
raw <- oligo::read.celfiles(huex_files)

# frma(raw, target = "core")  # Gene-level
# Output: ~22,000 rows (one per gene/transcript cluster)
# You get ONE value for entire TMEM129 gene

# Exon-level
exon_eset_frma <- frma(raw, target = "probeset")  # exon-level fRMA
exprs_exon <- Biobase::exprs(exon_eset_frma) # ; returns log2 scale
head(exprs_exon)


library(AnnotationDbi)
library(AnnotationHub)
library(huex10stprobeset.db)   # HuEx exon probeset → coords/symbol
library(ensembldb)
library(GenomicRanges); library(IRanges); library(S4Vectors)
library(AnnotationFilter)

ah <- AnnotationHub()
query_99 <- query(ah, c("EnsDb", "Homo sapiens", "Ensembl 99"))
query_99
edb <- query_99[[1]]

# exon coord
gr_exons <- exons(edb, return.type = "GRanges")
gr_exons

tmem_exons <- gr



# Build genomic ranges 

probe_ids <- rownames(exprs_exon)


# Annotate all probes
affy_ann <- AnnotationDbi::select(
  huex10stprobeset.db,
  keys    = probe_ids,
  keytype = "PROBEID",
  columns = c("SYMBOL","ENSEMBL","GENENAME")
)

# CRITICAL: Filter to unique mappings only
affy_ann_unique <- affy_ann %>%
  group_by(PROBEID) %>%
  dplyr::filter(n() == 1) %>%
  ungroup()

cat("Filtered from", nrow(affy_ann), "to", nrow(affy_ann_unique), "unique mappings\n")

# ENSG for TMEM129 you already derived as `ens_tmem`
tmem_probes <- subset(affy_ann, ENSEMBL %in% "ENSG00000168936")

# Subset your exon-level fRMA matrix to TMEM129 probesets
tmem_expr <- exprs_exon[intersect(tmem_probes$PROBEID, rownames(exprs_exon)), , drop = FALSE]
tmem_expr <-  as.data.frame(tmem_expr)
#tmem_expr$PROBEID <- rownames(tmem_expr)
# Merge info

tmem_probes <- tmem_probes %>% dplyr::select(PROBEID,SYMBOL)

tmem_expr <- tmem_expr %>%
  mutate(PROBEID = rownames(.)) %>%
  relocate(PROBEID, .before = 1)

final <- left_join(tmem_probes,tmem_expr,by="PROBEID")


# Note this pack coord is hg19
library(HuExExonProbesetLocation)

# Step 1: Get probeset genomic coordinates
# Simple filter - this is all you need!
probeset_coords <- HuExExonProbesetLocation %>%
  dplyr::filter(EPROBESETID %in% tmem_probes$PROBEID)

names(probeset_coords)[1] <- names(final)[1]

final_result <- left_join(probeset_coords,final,  by = "PROBEID")

```



