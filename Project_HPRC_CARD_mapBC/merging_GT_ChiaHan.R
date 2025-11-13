library(dplyr)
library(xlsx)



setwd('/Volumes/data/GWAS_COVNET_related/HPRC_COVNET_get_bams/Genotype/')

lsf <- list.files('.',pattern = '.txt')


# load the manifest
man <- read.csv('../../HPRC_ancestry_manifest.csv')
# make xx_1,xx_2 for chromosome
man$read_group <- sub("_[^_]+$", "", man$read_name)
man <- man[,c(-4)]

for (j in lsf){
  # j <- lsf[10]
  df <- read.delim(j,header =F,col.names = c('Sample',"chr","SNP_POS","GT"))
  df$rsid_with_chr <- paste0(gsub(".*(rs\\d+).*", "\\1", j),'_',df$chr,':',df$SNP_POS)
  
  df$Sample <- gsub('#','_',df$Sample)
  df$ID_person <- gsub('_.+','',df$Sample)
  df$read_group <- sub("_[^_]+$", "", df$Sample)
  names(df)[4] <- df$rsid_with_chr %>% unique()
  df.sub <- df[,c(6,7,4)]
  
  man <- man %>% left_join(df.sub,by=c('ID_person','read_group'))
  rm(df.sub)
  
}

#### Future uptades:####################################################
# some variants has DEL(*), need to replace to DEL
# double check some NAs. maybe realign with smaller 250 or 500K region or?
########################################################################

# SAVE
#write.table(man,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/HPRC and CARD genes regions/Genotyping from assemblies/HPRC_gt_result.csv',col.names = T,row.names = F,quote = F,sep = ',')


