# Goal: Counting Splicing Junction Reads for FGFR3 Exon Skipping and Normalization Methods

## Description:

### Reference Method from CHRNA5 Smoking Paper:

The quantification of splice junctions for all CHRNA5 exon 5 isoforms was performed using FeatureCounts (v4.5.0, in R) and a custom .bed file, which provided coordinates for each isoform-specific junction. Only samples with ≥30 total exon 5 splice junction reads were retained for analysis (A549: n=3, UMUC3: n=2; 67–211 junctions per sample); this excluded one UMUC3 replicate with insufficient junction read counts. Isoform ratios were calculated as the fraction of specific junctions over total junction reads. For each cell line, differences in isoform ratios between LT-CSCs and LT-DMSO were compared by unpaired two-sided t tests; p < 0.05 was considered significant.

**FeatureCounts command used by Max:**
```bash
##### The pET01_SJ.saf format: ##########
# Tab-delimited file:
# GeneID	Chr	Start	End	Strand
pET5p	pET01_CHRNA5_Exon5	1290	3593	+
##############################################
# featureCounts v2.0.6 
featureCounts -a pET01_SJ.saf -F SAF -o output input.bam --splitOnly -J
# The count output will have multiple columns that need to be annotated based on exon location
```

**Reference data:**
- Excel file location: `/LTG/Prokunina Lab/Paper Files/2025 CHRNA5 splicing/Submission 111025_Human Genomics/Hogshead Supplementary Tables_HG.xlsx`
- See **Table S4** for rationale and example workflow

### Additional Method: SH-SY5Y Junction Counting in CHRNA5 Region

```bash
# Using FeatureCounts to count junction reads (similar to Sashimi plot in IGV)
# System: ccad2
# Input: SH-SY5Y RNA-seq BAMs
# SAF file format (tab-separated):
# CHRNA5_iso	chr15	78588357	78593275	+

ml slurm/
ml subread/
for i in *.bam
do 
  featureCounts -t 12 -a Featurecount_CHRNA5_isoform/C5_iso.saf \
    -o Featurecount_CHRNA5_isoform/${i%.hg38.md.bam}_fct ${i} \
    --splitOnly -J -F SAF -p --countReadPairs
done
```

**Analysis code:**  
GitHub repository: [SHSY5Y_featurecount_CHRNA5.R](https://github.com/ChiaHan-1215/SH-SY5Y_cell_line_RNA_seq_with_long_read/blob/main/SHSY5Y_featurecount_CHRNA5.R)
