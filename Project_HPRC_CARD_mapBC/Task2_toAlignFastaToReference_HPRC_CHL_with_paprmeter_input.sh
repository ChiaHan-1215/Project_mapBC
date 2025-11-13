#!/bin/bash
#
# this file is myjob.sh
#
#SBATCH --job-name Task2-toAlignFastaToReference
#SBATCH --mail-type BEGIN,END
#
# load minimap2/2.26 and samtools/1.17
#module load minimap2/2.26
#module load samtools/1.17


## Author: Oscar Florez-Vargas


##########################
## == Set Parameters == ##

while getopts "I:h" opt; do
	case $opt in
		I) foldername="$OPTARG" ;;
		h) echo "  -I: folder full name"
		   exit 0 ;;
         \?) echo "Invalid option -$OPTARG" >&2
		     exit 1 ;;
	esac
done


genomicRegion=${foldername}

## Set path to where BAM files are being stored
PATH_TO_FILES=/data/leec20/GWAS_COVNET_related/HPRC_COVNET_get_bams/Retrive_regions
PATH_TO_REFERENCE=/data/Prokunina_Group/GenomeAssemblies/referenceData

##########################


# Create the "bam_files" folder
if [ ! -d $PATH_TO_FILES/$genomicRegion/bam_files ]; then
  mkdir -p $PATH_TO_FILES/$genomicRegion/bam_files
  chmod 775 $PATH_TO_FILES/$genomicRegion/bam_files;
fi


# Step 1: Concatenate all FASTA files in into input.fa
#for file in "$PATH_TO_FILES/$genomicRegion/fasta_files"/*.fa; do
#    filename=$(basename "$file" .fa)
#    sed -i "s/^>.*/>$filename/" "$file" > "$file.tmp" && mv "$file.tmp" "$file"
#done

cat $PATH_TO_FILES/$genomicRegion/fasta_files/*.fa > $PATH_TO_FILES/$genomicRegion/bam_files/HPRC_input.fa


# Step 2: Align the FASTA sequences to reference human genome build hg38
minimap2 \
	-Y -t 12 -a \
	$PATH_TO_REFERENCE/genome_hg38.fa \
	$PATH_TO_FILES/$genomicRegion/bam_files/HPRC_input.fa | \
	samtools sort -@ 6 -o $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam

# Step 3: Index the BAM file
samtools index $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam

# Step 4: Clean BAM file with only primary alignments
samtools view -F 0x900 -b $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.bam > $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.primary.bam
samtools index $PATH_TO_FILES/$genomicRegion/bam_files/HPRC.cb.primary.bam

# Step 5: remove the temporal FASTA file
# rm $OUTDIR/$genomicRegion/bam_files/HPRC_input.fa

