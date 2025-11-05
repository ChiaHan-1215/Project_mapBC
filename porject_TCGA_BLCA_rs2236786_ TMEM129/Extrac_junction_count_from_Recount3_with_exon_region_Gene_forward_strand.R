# Goal: extract junction count from Recount3 dataset. Target: TACC3 and FGFR3 
# Note: For this code is for gene that start from 5' to 3' -> direction. The gene from 3' to 5' is in other script 

library(recount3)
library(snapcount)
library(megadepth)
library(dplyr)
library(tidyr)

# the temp file is located in leec20/.cache/R/recount3
# needd to clean if too many?
# Note the TCGA file need to used in BIouwlf Rstudio, it's too large

human_projects <- available_projects()

# get the TCGA and GTEx 

#gtex <- subset(human_projects, file_source == "gtex" & project_type == "data_sources")
TCGA <- subset(human_projects, file_source == "tcga" & project_type == "data_sources")

list_TCGA <- TCGA$project


# for ( i in list_TCGA){

df <- data.frame()

i <- "BLCA"
target_tissue  <- TCGA %>% filter(project == i)

rse_jxn_b <- create_rse(
  target_tissue,
  type = "jxn"
)

# To find the TCGA ID, it's in the package 
TCGAIDmatch <- rse_jxn_b@colData@listData %>% as.data.frame()
TCGAIDmatch <- TCGAIDmatch[,c(2,4)]

matrix_data <- as.matrix(assay(rse_jxn_b))

df <- as.data.frame(matrix_data)
colnames(df) <- TCGAIDmatch$tcga.tcga_barcode[match(colnames(df), TCGAIDmatch$external_id)]



# TMEM129 region in hg38
# chr4:1715712-1721408
# chr4:1715936-1721342

# FGFR3 region in hg38
# chr4:1792585-1809586


df$POS <- rownames(df)

df_chr4 <- df[which(grepl("chr4:",df$POS)),]
df_chr4 <- df_chr4 %>%  separate(POS, into = c("chr", "start", "end",'strand'), sep = "[:-]")
df_chr4$strand <- gsub('chr4:[0-9]+-[0-9]+:','',rownames(df_chr4))

df_chr4$start <- as.numeric(df_chr4$start)
df_chr4$end <- as.numeric(df_chr4$end)

# df_chr4 <- df_chr4 %>% 
#   filter(start >= 1792585 & start <= 1809586 & 
#            end >= 1792585 & end <= 1809586)


df_chr4 <- df_chr4[,c('chr','start','end','strand',grep("TCGA",names(df_chr4),value = T))]
rownames(df_chr4) <- NULL


df_chr4$location <- paste0('chr4:',df_chr4$start,"-",df_chr4$end)
df_chr4 <- df_chr4[,c(grep("location",names(df_chr4),value = T),grep("TCGA",names(df_chr4),value = T))]

####### Get exons #######

library(rtracklayer)
library(dplyr)

gtf <- import('gencode.v39.annotation.gtf.gz')     

# c("FGFR3")

# 5 to 3 ->
genes_of_interest <-c("TACC3","FGFR3")
# 3 to 5 <-
# genes_of_interest <-c("FAM53A","SLBP","TMEM129")


# Since there's difference of 5' to 3' gene and 3' to 5' gene 
# so need to seperated by gene direct, the +1 -1 base for this code is for FGFR3,TACC3

# Subset exons for your genes
exons <- subset(gtf, type == "exon" & gene_name %in% genes_of_interest)
#exon_df <- as.data.frame(exons)[, c("seqnames", "start", "end", "gene_name", "exon_id", "exon_number")]
exon_df <- as.data.frame(exons)
# just get protein coding exon only
#exon_df <- exon_df[!is.na(exon_df$ccdsid),]
exon_df <- exon_df[, c("seqnames", "start", "end", "gene_name","transcript_name", "exon_id", "exon_number")]
exon_df$location <- paste0(exon_df$seqnames,':',exon_df$start,"-",exon_df$end)

# exon_df.u <- exon_df[!duplicated(exon_df[, c("exon_id", "location")]), ]

iso_list <- exon_df$transcript_name %>% unique() %>% sort()

jct_per_iso <- data.frame()


for (k in iso_list){
  
  # k <- iso_list[2]
  exon_df.u <- exon_df %>% filter(transcript_name == k)
  
  jct_df <- exon_df.u %>%
    arrange(transcript_name, as.numeric(exon_number)) %>%
    group_by(gene_name, transcript_name, seqnames) %>%
    mutate(
      next_start = lead(start),
      next_exon  = lead(exon_number)
    ) %>%
    filter(!is.na(next_start), as.numeric(next_exon) == as.numeric(exon_number) + 1) %>%
    ungroup() %>%    
    transmute(
      gene_name,
      transcript_name,
      jct_label = paste0("jct_ex", exon_number, "_ex", next_exon),
      location = paste0(seqnames, ":", end + 1L,"-", next_start - 1L)
    ) 
  jct_per_iso <- rbind(jct_per_iso,jct_df)
  
  
  
  
}


jct_per_iso <- jct_per_iso[!duplicated(jct_per_iso$location),]

df_final <- left_join(jct_per_iso,df_chr4,by='location')
# remove NA
df_final <- na.omit(df_final)

# Now save seperated by gene
df_final$gene_name %>% unique()

df_final.s <- df_final %>% filter(gene_name == "TACC3")

df_final.s <- df_final.s[,-c(1:3)]
df_final.s <- df_final.s[order(df_final.s$location),]
df_final.s.t <- as.data.frame(t(df_final.s))
names(df_final.s.t) <- df_final.s.t[1,]
df_final.s.t <- df_final.s.t[-1,]
df_final.s.t[] <- lapply(df_final.s.t[], as.numeric)
df_final.s.t <- df_final.s.t %>% mutate(TCGA_ID=rownames(df_final.s.t),.before = 1)

write.table(df_final.s.t,paste0("Recount_result/TCGA_BLCA_TACC3_jct_from_Recount3.csv"),col.names = T,row.names = F,sep = ',',quote = F)

rm(df_final.s)
rm(df_final.s.t)
