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

            
