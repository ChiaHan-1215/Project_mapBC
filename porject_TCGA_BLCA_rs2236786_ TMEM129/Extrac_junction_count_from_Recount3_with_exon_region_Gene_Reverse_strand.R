# Goal: extract junction count from Recount3 dataset. Target: Reverse strand gene direction (3 to 5) 

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



# TMEM129 - region in hg38
# chr4:1715712-1721408
# chr4:1715936-1721342

# FGFR3 + region in hg38
# chr4:1792585-1809586

# FAM53A -
# chr4:1636488-1687769

df$POS <- rownames(df)

df_chr4 <- df[which(grepl("chr4:",df$POS)),]
df_chr4 <- df_chr4 %>%  separate(POS, into = c("chr", "start", "end",'strand'), sep = "[:-]")
df_chr4$strand <- gsub('chr4:[0-9]+-[0-9]+:','',rownames(df_chr4))

df_chr4$start <- as.numeric(df_chr4$start)
df_chr4$end <- as.numeric(df_chr4$end)

# df_chr4 <- df_chr4 %>% 
#   filter(start >= 1636488 & start <= 1687769 & 
#            end >= 1636488 & end <= 1687769)


df_chr4 <- df_chr4[,c('chr','start','end','strand',grep("TCGA",names(df_chr4),value = T))]
rownames(df_chr4) <- NULL


df_chr4$location <- paste0('chr4:',df_chr4$start,"-",df_chr4$end)
df_chr4 <- df_chr4[,c(grep("location",names(df_chr4),value = T),grep("TCGA",names(df_chr4),value = T))]

####### Get exons #######

library(rtracklayer)
library(dplyr)

gtf <- import('gencode.v39.annotation.gtf.gz')     

# 5 to 3 ->
# genes_of_interest <-c("TACC3","FGFR3")
# 3 to 5 <-
genes_of_interest <-c("FAM53A","SLBP","TMEM129")


# Since there's difference of 5' to 3' gene and 3' to 5' gene 
# so need to seperated by gene direct, the +1 -1 base for this code is for FGFR3,TACC3

# Subset exons for your genes
exons <- subset(gtf, type == "exon" & gene_name %in% genes_of_interest)

# Make sure strand is included
exon_df <- as.data.frame(exons)[, c("seqnames","start","end","strand",
                                    "gene_name","transcript_name","exon_id","exon_number")]
exon_df$location <- paste0(exon_df$seqnames, ":", exon_df$start, "-", exon_df$end)

iso_list <- exon_df$transcript_name %>% unique() %>% sort()

jct_per_iso <- data.frame()

for (k in iso_list) {
  
  exon_df.u <- exon_df %>% filter(transcript_name == k)
  
  jct_df <- exon_df.u %>%
    arrange(as.numeric(exon_number)) %>%                             # transcript order
    group_by(gene_name, transcript_name, seqnames, strand) %>%
    mutate(
      next_start = lead(start),
      next_end   = lead(end),
      next_exon  = lead(exon_number)
    ) %>%
    filter(!is.na(next_start), as.numeric(next_exon) == as.numeric(exon_number) + 1) %>%
    mutate(
      # strand-aware junction bounds
      j_start = if_else(strand == "+", end + 1L,      next_end + 1L),
      j_end   = if_else(strand == "+", next_start - 1L, start   - 1L)
    ) %>%
    filter(j_start <= j_end) %>%                                      # drop overlaps/adjacent
    ungroup() %>%
    transmute(
      gene_name,
      transcript_name,
      jct_label = paste0("jct_ex", exon_number, "_ex", next_exon),
      location  = paste0(seqnames, ":", j_start, "-", j_end)
    )
  
  jct_per_iso <- rbind(jct_per_iso, jct_df)
}



jct_per_iso <- jct_per_iso[!duplicated(jct_per_iso$location),]

df_final <- left_join(jct_per_iso,df_chr4,by='location')
# remove NA
df_final <- na.omit(df_final)




# Now save seperated by gene
gs <- df_final$gene_name %>% unique()

df_final.s <- df_final %>% filter(gene_name == gs[3])

df_final.s <- df_final.s[,-c(1:3)]
df_final.s <- df_final.s[order(df_final.s$location),]
df_final.s.t <- as.data.frame(t(df_final.s))
names(df_final.s.t) <- df_final.s.t[1,]
df_final.s.t <- df_final.s.t[-1,]
df_final.s.t[] <- lapply(df_final.s.t[], as.numeric)
df_final.s.t <- df_final.s.t %>% mutate(TCGA_ID=rownames(df_final.s.t),.before = 1)

write.table(df_final.s.t,paste0("Recount_result/TCGA_BLCA_",gs[3],"_jct_from_Recount3.csv"),col.names = T,row.names = F,sep = ',',quote = F)

rm(df_final.s)
rm(df_final.s.t)
