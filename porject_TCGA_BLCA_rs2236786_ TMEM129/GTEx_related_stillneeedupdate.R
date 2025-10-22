
library(tidyverse)
library(dplyr)
library(ggsci)
library(ggrepel)
library(scales)
library(cowplot)
library(ggpubr)
library(rstatix)
library(patchwork)
library(tibble)
library(ggridges)
library(readr)


# load jct/transcript count from v10
# load Tissue mtadata
# load GT of rs223 

# mkae table of ID, full_ID, tissue, GT, jct, etc


df <- read.delim('~/Desktop/chr15_CHRNA5_bladder_project/UCSC_xena_TCGA_GTEx_dataset/CHRNA5_isosubset_isopct.txt')

# Now is to combined this to my GT and plot it 
# extract GTEx
df <- data.frame(t(df))
names(df) <- df[1,]
df$ID <- rownames(df)
df <- df[-1,]
rownames(df) <- NULL
df <- df[,c(3,2,1)]

df$ID <- gsub('\\.','-',df$ID)

df_GTEx <- df[which(grepl('GTEX',df$ID)),]
names(df_GTEx)[1] <- 'GTEx_ID'
# load the Tissue 

# newver <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/Metadata_Files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt')
setwd('~/Desktop/dbGAP_GTEX_and_other_dataset/GTEX_v9/TPM_and_rawcount_stuff/')
samid <- read.csv('../All_GTEx_sample_info_sex_age_race_etc/All_sample_info_donor_detailed/GTEx_Analysis_2021-02-11_v9_Annotations_GTEx_Analysis_2021-02-11_v9_Annotations_SampleAttributesDS.csv')
samid <- samid[,c(1,14)]
names(samid) <- c('GTEx_ID','tissue')




# merge, check duplicates

df_GTEx <- inner_join(df_GTEx,samid,by="GTEx_ID")
df_GTEx$PID <- gsub("(GTEX-[^-]*)-.*", "\\1", df_GTEx$GTEx_ID)

df_GTEx <- df_GTEx[,c(5,4,2,3)]


df_GTEx_test <- df_GTEx %>%
  group_by(PID, tissue) %>%
  summarise(
    ENST00000559554.5 = max(ENST00000559554.5),
    ENST00000299565.9 = max(ENST00000299565.9)
  )

df.wide <- df_GTEx_test %>% tidyr::pivot_wider(names_from = tissue, values_from = c(ENST00000559554.5,ENST00000299565.9))

# replace NA as No_data
#df.wide[] <- lapply(df.wide, function(x) replace(x, is.na(x), "No_data"))


# load Finaled data
ff <- read.csv('~/Desktop/chr15_CHRNA5_bladder_project/Finalized_GTEx_GT_isofrom_CHRAN5.csv')

# Merge with GT 
names(df.wide)[1] <- 'GTEx_ID'
ff <- ff[,c(1:83)]

final_TPM_perc <- inner_join(df.wide,ff,by="GTEx_ID")  
final_TPM_perc <- final_TPM_perc[,c(1,c(108:ncol(final_TPM_perc),2:107))]


final_TPM_perc$SEX <- as.factor(final_TPM_perc$SEX)
final_TPM_perc$RACE <- as.factor(final_TPM_perc$RACE)
final_TPM_perc$Smoking_status <- as.factor(final_TPM_perc$Smoking_status)

final_TPM_perc[, 84:ncol(final_TPM_perc)] <- lapply(final_TPM_perc[, 84:ncol(final_TPM_perc)], as.numeric)
final_TPM_perc[, which(grepl("_add",names(final_TPM_perc))) ] <- lapply(final_TPM_perc[, which(grepl("_add",names(final_TPM_perc)))], as.numeric)
final_TPM_perc <- data.frame(final_TPM_perc)
names(final_TPM_perc)[84:ncol(final_TPM_perc)] <- gsub('\\.\\.\\.|\\.\\.|\\.','_',names(final_TPM_perc)[84:ncol(final_TPM_perc)])
names(final_TPM_perc)[84:ncol(final_TPM_perc)] <- gsub('_$','',names(final_TPM_perc)[84:ncol(final_TPM_perc)])




df_list <- list()
df_list[["all_data"]] <- final_TPM_perc
df_list[["no_smoker"]] <- final_TPM_perc %>% filter(final_TPM_perc$Smoking_status == "No")
df_list[["yes_smoker"]] <- final_TPM_perc %>% filter(final_TPM_perc$Smoking_status == "Yes")
df_list[["smoke_over_15yrs"]] <- final_TPM_perc %>% filter(final_TPM_perc$Smoke_years > 15)
df_list[["smoke_less_15yrs"]] <- final_TPM_perc %>% filter(final_TPM_perc$Smoke_years <= 15 & final_TPM_perc$Smoking_status == 'Yes')
# df_list[['sk_muscle']] <- final_TPM_perc[,c(1:83,which(grepl("Muscle_Skeletal",names(final_TPM_perc))))]
# df_list[['sk_muscle_Yes_smoke']] <- final_TPM_perc[,c(1:83,which(grepl("Muscle_Skeletal",names(final_TPM_perc))))] %>% filter(final_TPM_perc$Smoking_status == 'Yes')
# df_list[['sk_muscle_No_smoke']] <- final_TPM_perc[,c(1:83,which(grepl("Muscle_Skeletal",names(final_TPM_perc))))] %>% filter(final_TPM_perc$Smoking_status == "No")


####
####
####

# ploting 


tissue <- gsub('ENST.+_[1-9]_','',names(final_TPM_perc)[84:ncol(final_TPM_perc)]) %>% unique()

# pdf('~/Desktop/chr15_CHRNA5_bladder_project/',height = 8,width = 10)

# df_sel <- df_list$no_smoker
# sk_color <- "#5FB2E8" # NO
# 
# df_sel <- df_list$yes_smoker
# sk_color <- "#FFBD95" # YES
# 
# df.out <- data.frame()
# df_count.summary <- data.frame()

# grep("_add",names(df_final),value = T)

# grep SNP 
for (i in c("rs142774214_add","rs503464_add","rs55853698_add","rs55781567_add")){
  # i <- "rs142774214_add"
  # i <- "rs503464_add"
  # i <- "rs55853698_add"
  # i <- "rs55781567_add"
  
  pdf(paste0('~/Desktop/chr15_CHRNA5_bladder_project/',i,'_Xena_isopct.pdf'),height = 9,width = 10)
  
  
  # grep tissue 
  for (k in tissue){
    # k <- tissue[2]
    df.out <- data.frame()
    df_count.summary <- data.frame()
    
    # dataset set as 
    df_sel <- df_list$no_smoker
    sk_color <- "#5FB2E8" # NO
    sm_st <- "No_Smoke"
    
    
    # k <- tissue[3]
    selected_col <- grep(k ,names(df_sel),value = T)
    tar <- df_sel[,c(names(df_sel)[1:7],gsub('_add',"",i),i,selected_col)]
    tar <- na.omit(tar)
    
    
    tmp.snp <- tar
    tmp.snp.0a <- tmp.snp[ which(tmp.snp[,i] == 0), ]
    tmp.snp.0b <- as.character(unique(tmp.snp.0a[gsub("_add", "", i)]))
    tmp.snp.1a <- tmp.snp[ which(tmp.snp[,i] == 1), ]
    tmp.snp.1b <- as.character(unique(tmp.snp.1a[gsub("_add", "", i)]))
    tmp.snp.2a <- tmp.snp[ which(tmp.snp[,i] == 2), ]
    tmp.snp.2b <- as.character(unique(tmp.snp.2a[gsub("_add", "", i)]))
    
    gender_sp <- c("Ovary" , 'Uterus' , 'Testis' ,  'Vagina',
                   "Prostate" , "Cervix_Ectocervix","Fallopian_Tube","Cervix_Endocervix")
    
    
    for (j in grep("ENST" ,names(tar),value = T)){
      if(any(gender_sp %in% k)){
        
        fmla <- as.formula(paste(j, "~", i, " + AGE"))
        
      } else { 
        
        fmla <- as.formula(paste(j, "~", i, " + AGE + SEX"))}
      
      out <- tryCatch(lm(fmla,tar), error=function(e){
        print(paste0('error of tissue ',k))
        return(NA)})
      
      summary(out)
      # 
      tmp.df <- data.frame(
        snp = i,
        variable = j,
        geno_0 = tmp.snp.0b,
        geno_1 = tmp.snp.1b,
        geno_2 = tmp.snp.2b,
        n_0 = nrow(tmp.snp.0a),
        n_1 = nrow(tmp.snp.1a),
        n_2 = nrow(tmp.snp.2a),
        
        beta = coef(summary(out))[2,1],
        p_val = coef(summary(out))[2,4])
      # tmp.df$variable <- paste0(k,"_",tmp.df$variable)
      
      rownames(tmp.df) <-NULL
      
      df.out <- rbind(df.out, tmp.df)
      
      ####
    }
    
    
    FLbp <- df.out %>% filter(variable == paste0("ENST00000299565_9_",k))
    SSbp <- df.out %>% filter(variable == paste0("ENST00000559554_5_",k))
    
    
    
    p1 <- ggplot(tar, aes(x=as.factor(.data[[i]]), y=.data[[paste0("ENST00000299565_9_",k)]])) +
      geom_violin(width = 0.5, lwd = 0.05, colour = "black",fill=sk_color,
                  show.legend = NA, inherit.aes = TRUE) +
      geom_boxplot(width = 0.04, lwd = 0.3, colour = "black",
                   outlier.colour = NA, outlier.shape = NA,
                   show.legend = NA, inherit.aes = TRUE, fill = "white") +
      geom_point(
        size=0.5,
        alpha=0.3, position = position_jitter(width=0.05,height = 0.01,seed = 100))+
      
      stat_summary(fun=median, geom="point", shape=20, size=2, color="#FF0000") +
      theme_classic() +
      theme(legend.position="none",axis.title.x=element_blank()) + #xlab(paste0(i,' , ',k)) +
      annotate("text", x = Inf, y = Inf, label = paste0(sm_st," beta: ",round(FLbp[9],3),", p-val: ", round(FLbp[10],3)), hjust = 1, vjust = 1)
    
    
    
    #### print SS
    
    p2 <- ggplot(tar, aes(x=as.factor(.data[[i]]), y=.data[[paste0("ENST00000559554_5_",k)]])) +
      geom_violin(width = 0.5, lwd = 0.05, colour = "black",fill=sk_color,
                  show.legend = NA, inherit.aes = TRUE) +
      geom_boxplot(width = 0.04, lwd = 0.3, colour = "black",
                   outlier.colour = NA, outlier.shape = NA,
                   show.legend = NA, inherit.aes = TRUE, fill = "white") +
      geom_point(
        size=0.5,
        alpha=0.3, position = position_jitter(width=0.05,height = 0.01,seed = 100))+
      
      stat_summary(fun=median, geom="point", shape=20, size=2, color="#FF0000") +
      theme_classic() +
      theme(legend.position="none",axis.title.x=element_blank()) + #xlab(paste0(i,',',k)) +
      annotate("text", x = Inf, y = Inf, label = paste0(sm_st," beta: ",round(SSbp[9],3),", p-val: ", round(SSbp[10],3)), hjust = 1, vjust = 1)
    
    
    
    
    ###### ###### ###### ######
    ###### ###### ###### ######
    ##### SMOKING data ######
    ###### ###### ###### ######
    ###### ###### ###### ######
    
    df_sel <- df_list$yes_smoker
    sk_color <- "#FFBD95" # YES
    sm_st <- "Yes_Smoke"
    
    df.out <- data.frame()
    df_count.summary <- data.frame()
    
    selected_col <- grep(k ,names(df_sel),value = T)
    tar <- df_sel[,c(names(df_sel)[1:7],gsub('_add',"",i),i,selected_col)]
    tar <- na.omit(tar)
    
    
    tmp.snp <- tar
    tmp.snp.0a <- tmp.snp[ which(tmp.snp[,i] == 0), ]
    tmp.snp.0b <- as.character(unique(tmp.snp.0a[gsub("_add", "", i)]))
    tmp.snp.1a <- tmp.snp[ which(tmp.snp[,i] == 1), ]
    tmp.snp.1b <- as.character(unique(tmp.snp.1a[gsub("_add", "", i)]))
    tmp.snp.2a <- tmp.snp[ which(tmp.snp[,i] == 2), ]
    tmp.snp.2b <- as.character(unique(tmp.snp.2a[gsub("_add", "", i)]))
    
    gender_sp <- c("Ovary" , 'Uterus' , 'Testis' ,  'Vagina',
                   "Prostate" , "Cervix_Ectocervix","Fallopian_Tube","Cervix_Endocervix")
    
    
    for (j in grep("ENST" ,names(tar),value = T)){
      if(any(gender_sp %in% k)){
        
        fmla <- as.formula(paste(j, "~", i, " + AGE"))
        
      } else { 
        
        fmla <- as.formula(paste(j, "~", i, " + AGE + SEX"))}
      
      out <- tryCatch(lm(fmla,tar), error=function(e){
        print(paste0('error of tissue ',k))
        return(NA)})
      
      summary(out)
      # 
      tmp.df <- data.frame(
        snp = i,
        variable = j,
        geno_0 = tmp.snp.0b,
        geno_1 = tmp.snp.1b,
        geno_2 = tmp.snp.2b,
        n_0 = nrow(tmp.snp.0a),
        n_1 = nrow(tmp.snp.1a),
        n_2 = nrow(tmp.snp.2a),
        
        beta = tryCatch(coef(summary(out))[2,1], error=function(e){
          print(paste0('error of tissue ',k))
          return(NA)}) ,
        p_val = tryCatch(coef(summary(out))[2,4], error=function(e){
          print(paste0('error of tissue ',k))
          return(NA)})) 
      
      
      rownames(tmp.df) <-NULL
      
      df.out <- rbind(df.out, tmp.df)
      
      ####
    }
    
    
    
    # ENST00000299565
    FLbp <- df.out %>% filter(variable == paste0("ENST00000299565_9_",k))
    SSbp <- df.out %>% filter(variable == paste0("ENST00000559554_5_",k))
    
    
    
    p3 <- ggplot(tar, aes(x=as.factor(.data[[i]]), y=.data[[paste0("ENST00000299565_9_",k)]])) +
      geom_violin(width = 0.5, lwd = 0.05, colour = "black",fill=sk_color,
                  show.legend = NA, inherit.aes = TRUE) +
      geom_boxplot(width = 0.04, lwd = 0.3, colour = "black",
                   outlier.colour = NA, outlier.shape = NA,
                   show.legend = NA, inherit.aes = TRUE, fill = "white") +
      geom_point(
        size=0.5,
        alpha=0.3, position = position_jitter(width=0.05,height = 0.01,seed = 100))+
      
      stat_summary(fun=median, geom="point", shape=20, size=2, color="#FF0000") +
      theme_classic() +
      theme(legend.position="none") + xlab(paste0(i,' , ',k)) +
      annotate("text", x = Inf, y = Inf, label = paste0(sm_st, " beta: ",round(FLbp[9],3),", p-val: ", round(FLbp[10],3)), hjust = 1, vjust = 1)
    
    
    
    #### print SS
    
    p4 <- ggplot(tar, aes(x=as.factor(.data[[i]]), y=.data[[paste0("ENST00000559554_5_",k)]])) +
      geom_violin(width = 0.5, lwd = 0.05, colour = "black",fill=sk_color,
                  show.legend = NA, inherit.aes = TRUE) +
      geom_boxplot(width = 0.04, lwd = 0.3, colour = "black",
                   outlier.colour = NA, outlier.shape = NA,
                   show.legend = NA, inherit.aes = TRUE, fill = "white") +
      geom_point(
        size=0.5,
        alpha=0.3, position = position_jitter(width=0.05,height = 0.01,seed = 100))+
      
      stat_summary(fun=median, geom="point", shape=20, size=2, color="#FF0000") +
      theme_classic() +
      theme(legend.position="none") + xlab(paste0(i,',',k)) +
      annotate("text", x = Inf, y = Inf, label = paste0(sm_st, " beta: ",round(SSbp[9],3),", p-val: ", round(SSbp[10],3)), hjust = 1, vjust = 1)
    
    
    p <- p1 + p2 + p3 + p4 # + ggtitle("SMOKING_NO")
    
    print(p)
    
    
  }
  
  dev.off()
  
}



