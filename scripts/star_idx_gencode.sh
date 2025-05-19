#!/bin/bash
#SBATCH --job-name=star_idx_gencode
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 17
#SBATCH --mem 128G
#SBATCH -t 8:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=james.ri
#SBATCH --out=/home/%u/logs/%x_%j.log
#SBATCH --error=/home/%u/logs/%x_%j.err

shopt -s expand_aliases
CRN=$( basename $( find /courses -maxdepth 1 -type d -name "BINF6430.*" ) | cut -d "." -f2 )

alias STAR='singularity run -B "/courses:/courses,/home:/home,/scratch:/scratch" /shared/container_repository/explorer/star/2.7.10b/star_2.7.10b.sif STAR'

BASE_DIR=/courses/BINF6430.${CRN}
DATA_DIR=${BASE_DIR}/shared/gencode

GENOME_DIR=${BASE_DIR}/students/${USER}/STAR_index  # Changed from STAR to STAR_index
mkdir -p ${GENOME_DIR}

# STAR genome indexing command
STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${DATA_DIR}/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile ${DATA_DIR}/gencode.v46.annotation.gtf \
     --sjdbOverhang 100

echo "STAR genome indexing completed."
