# Goal: retrive region of BAMs around 1Mb of variant of interest


### Step1
- The list is located in `/Volumes/Prokunina_Group/GenomeAssemblies/mapBC_GWAS_signals_09302025.xlsx`
- the SNPs is in hg19, but need to change to hg38, then expand to left/right 500kb so that total is 1mb


- **FOR HPRC**, The script to follow is `/data/Prokunina_Group/GenomeAssemblies/scripts_toRetrieve_HPRC_regions_fromBiowulf/Task1_RetrieveGenomeRegion_HPRCr2_COVNET_10p23.31.sh`
- Use script for guide, chagne input value like POS CHR SNP etc
- named the folder as mapBC_rsID_cytoband (ex: mapBC_rs1234_10p23.31 something like that)

- **FOR CARD**, The script to follow is `/data/Prokunina_Group/GenomeAssemblies/scripts_toRetrieve_CARD_regions_fromLocal/Task1_RetrieveGenomeRegion_CARD_ztwo_COVNET_10p23.31.sh`
- use same stragty above

*************

### Setp2

Step 2 is mainly alignment, the script to follow is `/data/Prokunina_Group/GenomeAssemblies/Task2_toAlignFastaToReference_CARD_COVNET_10p23.31.sh`

- modfy nesessary value and run it

*************

### Step3 

Calling variants, still under Oscar constructiuon
