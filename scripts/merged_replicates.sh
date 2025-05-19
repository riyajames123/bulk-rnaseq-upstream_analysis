#!/bin/bash
#SBATCH --job-name=merge_replicates
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 17
#SBATCH --mem 128G
#SBATCH -t 8:00:00
#SBATCH --array=1-12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=m.sherman@northeastern.edu
#SBATCH --out=/home/%u/logs/%x_%j_%A_%a.log
#SBATCH --error=/home/%u/logs/%x_%j_%A_%a.err
shopt -s expand_aliases
CRN=$( basename $( find /courses -maxdepth 1 -type d -name "BINF6430.*" ) | cut -d "." -f2 )

module load singularity/3.10.3
alias samtools="singularity run -B '/courses:/courses' /courses/BINF6430.${CRN}/shared/singularity_containers/samtools-latest.sif samtools"

BASE=/courses/BINF6430.${CRN}
DATA=${BASE}/data/ReaganData/aligned_gencode

ACC_FILE=${BASE}/data/ReaganData/merged_replicates.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ACC_FILE})
SAMPLE_SUFFIX="Aligned.toTranscriptome.out.bam"
LEFT=${DATA}/${SAMPLE}_L001_${SAMPLE_SUFFIX}
RIGHT=${DATA}/${SAMPLE}_L002_${SAMPLE_SUFFIX}

SAMPLE_OUTPUT="${SAMPLE}_merged.aligned.toTranscriptome.out.bam"

OUT_DIR=${DATA}/merged
mkdir -p ${OUT_DIR}

samtools merge -@ ${SLURM_CPUS_PER_TASK} \
    ${OUT_DIR}/${SAMPLE_OUTPUT} \
    -fO BAM \
    ${LEFT} \
    ${RIGHT}