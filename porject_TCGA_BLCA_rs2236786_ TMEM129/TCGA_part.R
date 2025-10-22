#TCGA 
#code extracted from /Desktop/TERT_R2_project/TR2_manuscript_scripts_and_data/GTEX_TCGA_TPM_beta.FL_plot_v2.R


library(dplyr)
library(ggplot2)
library(ggridges)
library(hrbrthemes)
library(viridis)
library(ggrain)
library(ggthemes)
library(see)
library(patchwork)


# on L-drive
setwd('/Volumes/LTG/Prokunina Lab/Paper Files/2023 TERT_T2 paper/final datasets/TCGA_GTEx_isoform_plotting_01142025/Original_source_data_and_script/')


################################################################################
# TCGA part 
################################################################################


ta <- read.csv('tcga_TERT_rsem_isoform_tpm.txt',header = T,sep = '\t')
options(scipen = 999)
rownames(ta) <- ta$sample
ta <- ta[,-1]
ta <- ta %>%  mutate_all(~ 2^. - 0.001)
ta <- apply(ta, MARGIN = c(1, 2), FUN = function(x) round(x, 2))

ta <- data.frame(t(ta))
ta$total_TERT <- rowSums(ta[,c(1:7)])
ta <- ta[,c("ENST00000310581.9","ENST00000508104.2",'ENST00000484238.6','total_TERT')]
ta$PID <- rownames(ta)
ta$PID <- gsub('\\.','-',ta$PID)
# load sample info

infota <- read.table('TCGA_phenotype_denseDataOnlyDownload.tsv',header = T,sep = '\t')
names(infota)[1] <- 'PID'

infota$TCGA_cancer_id <- car::recode(infota$X_primary_disease," 'acute myeloid leukemia'='LAML' ; 'adrenocortical cancer'='ACC' ; 'bladder urothelial carcinoma'='BLCA' ; 'brain lower grade glioma'='LGG' ; 'breast invasive carcinoma'='BRCA' ; 'cervical & endocervical cancer'='CESC';
'cholangiocarcinoma'='CHOL' ; 'colon adenocarcinoma'='COAD' ; 'diffuse large B-cell lymphoma'='DLBC' ; 'esophageal carcinoma'='ESCA' ; 'glioblastoma multiforme'='GBM' ; 'head & neck squamous cell carcinoma'='HNSC' ;
'kidney chromophobe'='KICH' ; 'kidney clear cell carcinoma'='KIRC' ; 'kidney papillary cell carcinoma'='KIRP' ; 'liver hepatocellular carcinoma'='LIHC' ; 'lung adenocarcinoma'='LUAD' ;
'lung squamous cell carcinoma'='LUSC' ; 'mesothelioma'='MESO' ; 'ovarian serous cystadenocarcinoma'='OV' ; 'pancreatic adenocarcinoma'='PAAD' ; 'pheochromocytoma & paraganglioma'='PCPG' ;
'prostate adenocarcinoma'='PRAD' ; 'rectum adenocarcinoma'='READ' ; 'sarcoma'='SARC' ; 'skin cutaneous melanoma'='SKCM' ; 'stomach adenocarcinoma'='STAD' ;
'testicular germ cell tumor'='TGCT' ; 'thymoma'='THYM' ; 'thyroid carcinoma'='THCA' ; 'uterine carcinosarcoma'='UCS' ; 'uterine corpus endometrioid carcinoma'='UCEC' ;
                                     'uveal melanoma'='UVM' ")

# join

tacb <- inner_join(infota,ta,by = 'PID')

### the final is to filterout the Normal part!!!
tacb <- tacb[!grepl('Normal',tacb$sample_type),]



############### From Oscar's script ###############

library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(data.table)

setwd("/DCEG/Branches/LTG/Prokunina/TCGA_data")
getwd()

df.id <- read.table(file = "All_TCGA_case_IDs/TCGA_person_IDs_withCancerType_fromGDC.txt", sep = "\t", header = TRUE)
table(df.id$Cancer.Type)



## Bulk isoform expression, TPM in hg38

df_iso.transcript <- read.table(file = "Original source data Expression/UCSC Toil RNA-seq Recompute/transcript expression RNAseq_hg38/tcga_rsem_isoform_tpm.gz", header = TRUE, sep = "\t")
df_iso.name <- read.table(file = "Original source data Expression/UCSC Toil RNA-seq Recompute/transcript expression RNAseq_hg38/probeMap2Fgencode.v23.annotation.transcript.probemap", header = TRUE, sep = "\t")
            
## gene expression, TPM in hg38
df.gene <- read.table(file = "Original source data Expression/UCSC Toil RNA-seq Recompute/gene expression RNAseq_hg38/tcga_RSEM_gene_tpm.gz", header = TRUE, sep = "\t")
df_gene.name <- read.table(file = "Original source data Expression/UCSC Toil RNA-seq Recompute/gene expression RNAseq_hg38/probeMap2Fgencode.v23.annotation.gene.probemap", header = TRUE, sep = "\t")



            
names(df_iso.transcript) <- gsub("[.]", "-", names(df_iso.transcript))
head(df_iso.transcript[1:10, 1:5])
head(df_iso.name[1:10, 1:5])


# importing and formatting 
df.iso.exp <- merge(df_iso.name, df_iso.transcript, by.x = "id", by.y = "sample")
head(df.iso.exp[1:10, 1:9])

# formatting TCGA barcode
tcga_code <- as.data.frame(names(df_iso.transcript))
names(tcga_code)[1] <- "TCGA_barcode"
tcga_code$ID <- tcga_code$TCGA_barcode
tcga_code <- tcga_code %>% separate(ID, into = paste("ID", 1:7, sep = "_"))
tcga_code$sample <- as.numeric(str_extract(tcga_code$ID_4, "[0-9]+"))
tcga_code$ID_person <- do.call(paste, c(tcga_code[c("ID_1", "ID_2", "ID_3")], sep = "-"))
tcga_code$TCGA_barcode <- do.call(paste, c(tcga_code[c("ID_1", "ID_2", "ID_3", "ID_4")], sep = "-"))
tcga_code <- tcga_code[ , c(1,10,9)]
tcga_code <- merge(df.id, tcga_code, by = "ID_person", all = TRUE)
tcga_code$sample <- replace(tcga_code$sample, tcga_code$Cancer.Type != "SKCM" & tcga_code$sample == 1 | tcga_code$sample == 3 |
                              tcga_code$Cancer.Type == "SKCM" & tcga_code$sample == 6, "tumor")
tcga_code$sample <- replace(tcga_code$sample, tcga_code$sample == 11, "normal")

tcga_code <- tcga_code[ which(tcga_code$sample == "tumor" | tcga_code$sample == "normal"), ]
table(tcga_code$sample)
table(tcga_code$Cancer.Type)
tcga_code <- tcga_code[ complete.cases(tcga_code$Cancer.Type), ]


for ( i in unique(tcga_code$Cancer.Type) ) {
  # i <- "BLCA"
  
  tmp.ca <- df_iso.transcript[ , c(names(df_iso.transcript[1]), tcga_code[ which(tcga_code$Cancer.Type == i), ]$TCGA_barcode ) ]
  row.names(tmp.ca) <- tmp.ca$sample
  
  df.ca.exp <- as.data.frame(t(tmp.ca))
  head(df.ca.exp[1:10,1:10])
  
  colnames(df.ca.exp) <- as.character(unlist(df.ca.exp[1,]))
  df.ca.exp <- df.ca.exp[-1, ]
  df.ca.exp$TCGA_barcode <- rownames(df.ca.exp)
  df.ca.exp$TCGA_barcode <- gsub("[.].", "", df.ca.exp$TCGA_barcode)
  row.names(df.ca.exp) <- NULL
  n_occur <- as.data.frame(table(df.ca.exp$TCGA_barcode))
  df.ca.exp <- df.ca.exp[!duplicated(df.ca.exp$TCGA_barcode), ]
  df.ca.exp$TCGA_barcode
  
  row.names(df.ca.exp) <- df.ca.exp$TCGA_barcode
  
  # eliminate duplicated samples and get only tumor samples
  tcga_tumor <- tcga_code[!duplicated(tcga_code[c("ID_person", "sample")]) & tcga_code$sample == "tumor", ]
  tcga_tumor <- tcga_tumor[!is.na(tcga_tumor$ID_person),]
  table(tcga_tumor$Cancer.Type)
  
  # eliminate duplicated samples and get only normal samples
  tcga_normal <- tcga_code[!duplicated(tcga_code[c("ID_person", "sample")]) & tcga_code$sample == "normal", ]
  tcga_normal <- tcga_normal[!is.na(tcga_normal$ID_person),]
  table(tcga_normal$Cancer.Type)
  
  # get working data
  df.ISO.raw.tumor <- merge(tcga_tumor, df.ca.exp, by = "TCGA_barcode")
  df.ISO.raw.normal <- merge(tcga_normal, df.ca.exp, by = "TCGA_barcode")
  table(df.ISO.raw.tumor$Cancer.Type)
  table(df.ISO.raw.normal$Cancer.Type)
  
  # get subsets for tumor samples
  row.names(df.ISO.raw.tumor) <- df.ISO.raw.tumor$ID_person
  df.ISO.raw.tumor <- data.frame(t(df.ISO.raw.tumor[,5:ncol(df.ISO.raw.tumor)]))
  names(df.ISO.raw.tumor) <- gsub("[.]", "-", names(df.ISO.raw.tumor))
  df.ISO.raw.tumor <- merge(df_iso.name, df.ISO.raw.tumor, by.x = "id", by.y = "row.names")
  
  ## get subsets for normal samples
  # row.names(df.ISO.raw.normal) <- df.ISO.raw.normal$ID_person
  # df.ISO.raw.normal <- data.frame(t(df.ISO.raw.normal[,5:ncol(df.ISO.raw.normal)]))
  # names(df.ISO.raw.normal) <- gsub("[.]", "-", names(df.ISO.raw.normal))
  # df.ISO.raw.normal <- merge(df_iso.name, df.ISO.raw.normal, by.x = "id", by.y = "row.names")
  
  pathToFile = "Original source data Expression/UCSC Toil RNA-seq Recompute/transcript expression RNAseq_hg38/byCancerType"
  
  ## Saving data in TXT files for re-runnning from here
  write.table(df.ISO.raw.tumor, file = paste0(pathToFile, "/TCGA_", i, "_IsoformExpression.Profile_OnlyTumor.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  # write.table(df.ISO.raw.normal, file = paste0(pathToFile, "/TCGA_", i, "_IsoformExpression.Profile_OnlyNormal.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
}

            
