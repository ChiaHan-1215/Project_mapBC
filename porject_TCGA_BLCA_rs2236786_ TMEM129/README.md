## Goal: TMEM129 related 
- TCGA BLCA T/N rs2236786 GT in T-drive:
- Get junction readcount TMEM129 from TCGA data from Recount3 package

## Code: 

`Extrac_junction_count_from_Recount3_with_exon_region.R` -  How to extract junction count from target gene, need to run in Biowulf. Add more details of how to get each exons region in order to get intron region. Get all the isofroms junction region and align it to the Recount3 region


- **IMVIGOR210 dataset**

## make markdown or notes for how we download and get z-score/raw count etc from IMVagor210 dataset 




## To do list

- finding TPM of TCGA/GTEx of TMEM129 isoform, can we use UCSC Xena data we dowlaoded for TERT paper?

  The locatation of TCGA file:
`/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/Original source data Expression/UCSC Toil RNA-seq Recompute`

 - GTEx v10 have TPM,read count, junction count etc in `https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression`

  We also have folder for v10: 

  `/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10`



**check if we had most of file related to TCGA immune https://docs.google.com/document/d/16ikNbWpkIi_oTR3mYVu-Mle8fyLxxcgV/edit**

the link: https://gdc.cancer.gov/about-data/publications/panimmune

- if not, download file in T-drive and make note that file is downloaded

```
# In Biowulf

ml gdc-client
gdc-client download -t token.txt file_uuid/id

```
## Once dowloaded all MHC related file, see what's in the data and organizing into masterfile

 





