#!/bin/bash
#SBATCH --job-name=star_gencode_align
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 17
#SBATCH --mem=128G
#SBATCH -t 8:00:00
#SBATCH --array=1-12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=james.ri@northeastern.edu
#SBATCH --out=/home/%u/logs/%x_%j.log
#SBATCH --error=/home/%u/logs/%x_%j.err

shopt -s expand_aliases
CRN=$(basename $(find /courses -maxdepth 1 -type d -name "BINF6430.*") | cut -d "." -f2)

# Load STAR alias
alias STAR='singularity run -B "/courses:/courses,/home:/home,/scratch:/scratch" /shared/container_repository/explorer/star/2.7.10b/star_2.7.10b.sif STAR'

# Paths
BASE_DIR=/courses/BINF6430.${CRN}
USER_DIR=${BASE_DIR}/students/${USER}
GENOME_DIR=${USER_DIR}/STAR_index
GTF_FILE=${BASE_DIR}/shared/gencode/gencode.v46.annotation.gtf

OUTPUT_DIR=${USER_DIR}/data/reagan_aligned
mkdir -p ${OUTPUT_DIR}
TMP_DIR=${USER_DIR}/data/reagan_aligned_tmp
mkdir -p ${TMP_DIR}

# Sample Directory
SAMPLE_BASE=${BASE_DIR}/data/ReaganData/Reagan_PE85_TakaraPicoV2_HC_CM_10042022
SAMPLES=($(ls -d ${SAMPLE_BASE}/*/))

# Select sample based on SLURM_ARRAY_TASK_ID
SAMPLE_DIR=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
SAMPLE_NAME=$(basename ${SAMPLE_DIR})

echo "Processing sample: ${SAMPLE_NAME}"

# Identify Read Files (Merging Lanes)
LEFT_MERGED=${TMP_DIR}/${SAMPLE_NAME}_R1.fastq.gz
RIGHT_MERGED=${TMP_DIR}/${SAMPLE_NAME}_R2.fastq.gz

cat ${SAMPLE_DIR}/*_L001_R1_001.fastq.gz ${SAMPLE_DIR}/*_L002_R1_001.fastq.gz > ${LEFT_MERGED}
cat ${SAMPLE_DIR}/*_L001_R2_001.fastq.gz ${SAMPLE_DIR}/*_L002_R2_001.fastq.gz > ${RIGHT_MERGED}

# Check if files exist
if [[ ! -s "$LEFT_MERGED" || ! -s "$RIGHT_MERGED" ]]; then
    echo "Skipping $SAMPLE_NAME (No R1/R2 FASTQ files found after merging)"
    exit 1
fi

echo "Merged lanes for ${SAMPLE_NAME}, starting alignment..."

STAR --runThreadN 17 \
     --genomeDir ${GENOME_DIR} \
     --sjdbGTFfile ${GTF_FILE} \
     --readFilesIn ${LEFT_MERGED} ${RIGHT_MERGED} \
     --readFilesCommand zcat \
     --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE_NAME}_ \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 100000000000

echo "Alignment complete for ${SAMPLE_NAME}!"
