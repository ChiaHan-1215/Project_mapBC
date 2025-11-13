#!/bin/bash

#SBATCH --job-name RetrieveGenomicRegion
#SBATCH --mail-type BEGIN,END

# Authors: Tongwu Zhang & Oscar Florez-Vargas

module load seqkit

##########################
## == Set Parameters == ##

# Default values
runmod="local"
nthread=10
build="hg38"

# Parse command line arguments
#The format bc_snp_gene_cytoband 

while getopts "c:s:e:r:t:b:g:h" opt; do
  case $opt in
    c) chrname="$OPTARG" ;;
    s) Position1="$OPTARG" ;;
    e) Position2="$OPTARG" ;;
    r) runmod="$OPTARG" ;;
    t) nthread="$OPTARG" ;;
    b) build="$OPTARG" ;;
    g) taggene="$OPTARG" ;;
    h) echo "Usage: sbatch $0 -c <chr> -s <start_pos> -e <end_pos> -g <taggene> [-r <runmod>] [-t <nthread>] [-b <build>]"
       echo "  -c: Chromosome name (e.g., chr10)"
       echo "  -s: Start position"
       echo "  -e: End position"
       echo "  -g: Tag gene name"
       echo "  -r: Run mode (default: local)"
       echo "  -t: Number of threads (default: 10)"
       echo "  -b: Genome build (default: hg38)"
       exit 0 ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1 ;;
  esac
done



# Check required arguments
if [ -z "$chrname" ] || [ -z "$Position1" ] || [ -z "$Position2" ] || [ -z "$taggene" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch $0 -c <chr> -s <start_pos> -e <end_pos> -g <taggene> [-r <runmod>] [-t <nthread>] [-b <build>]"
    exit 1
fi


genomicRegion="${chrname}_${Position1}_${Position2}"
storageFolder="${build}_${chrname}_${Position1}_${Position2}_${taggene}"

reference_fasta="/data/Prokunina_Group/GenomeAssemblies/referenceData/genome_hg38.fa"
PATH_TO_FILES="/data/Prokunina_Group/GenomeAssemblies/CARD/assemblies"
OUTDIR="/data/Prokunina_Group/GenomeAssemblies/CARD/retrievedRegions"

# Create the storage folder
if [ ! -d $OUTDIR/$storageFolder ]; then
  mkdir -p $OUTDIR/$storageFolder
  chmod 775 $OUTDIR/$storageFolder;
fi

# Create the "fasta_files" folder
if [ ! -d $OUTDIR/$storageFolder/fasta_files ]; then
  mkdir -p $OUTDIR/$storageFolder/fasta_files
  chmod 775 $OUTDIR/$storageFolder/fasta_files;
fi

working_dir=$OUTDIR/$storageFolder/fasta_files

PrimerLen=25
ExtendLen=200
nmism=0

# Default to current directory if no working directory is provided
if [ -z "$working_dir" ]; then
    working_dir=$(pwd)
fi

cd $working_dir


# define External Primers
PrimerPos1=$((Position1 + PrimerLen))
PrimerPos2=$((Position2 - PrimerLen))

Primer_left=$(seqkit faidx $reference_fasta ${chrname}:${Position1}-${PrimerPos1} | tail -1)
Primer_right=$(seqkit faidx $reference_fasta ${chrname}:${PrimerPos2}-${Position2} | seqkit seq --seq-type DNA -r -p | tail -1)

# define External Primers
Position1_exd=$((Position1 - ExtendLen))
Position2_exd=$((Position2 + ExtendLen))
PrimerPos1_exd=$((Position1_exd + PrimerLen))
PrimerPos2_exd=$((Position2_exd - PrimerLen))

Primer_left_exd=$(seqkit faidx $reference_fasta ${chrname}:${Position1_exd}-${PrimerPos1_exd} | tail -1)
Primer_right_exd=$(seqkit faidx $reference_fasta ${chrname}:${PrimerPos2_exd}-${Position2_exd} | seqkit seq --seq-type DNA -r -p | tail -1)

# run in the local
if [ "$runmod" = "local" ]; then
  for i in $(ls "$PATH_TO_FILES"/*.fa.gz); do
    filename=$(basename "$i" .fa.gz)
    time seqkit amplicon -w 0 -m $nmism -M -j $nthread -F "$Primer_left" -R "$Primer_right" "$i" >"${filename}.fa"  
    if [ ! -s "${filename}.fa" ]; then
      time seqkit amplicon -w 0 -m $nmism -M -j 6 -F "$Primer_left_exd" -R "$Primer_right_exd" "$i" | seqkit subseq -r ${ExtendLen}:-${ExtendLen}  >"${filename}.fa"
    fi
  done
  
  ## delete the empty file
  find . -maxdepth 1 -type f -empty -delete

else
  ls "$PATH_TO_FILES"/*.fa.gz | /data/Prokunina_Group/GenomeAssemblies/software/rush --dry-run "seqkit amplicon -w 0 -m $nmism -M -j $nthread -F $Primer_left -R $Primer_right {} >{%..}.fa; if [ ! -s {%..}.fa ]; then seqkit amplicon -w 0 -m $nmism -M -j 6 -F $Primer_left_exd -R $Primer_right_exd {} | seqkit subseq -r ${ExtendLen}:-${ExtendLen}  >{%..}.fa; else true; fi" >cmdfile_$genomicRegion

  #submit job to biowulf
  swarm -f cmdfile_$genomicRegion -t $nthread -g 80 --logdir logs_$genomicRegion  --job-name $genomicRegion --time=3:00:00 --module seqkit --partition quick 
fi

