#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --partition=courses
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem 128G
#SBATCH -t 8:00:00
#SBATCH --array=1-12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=james.ri@northeastern.edu
#SBATCH --output=/courses/BINF6430.202530/students/james.ri/logs/salmon_quant_%j.out
#SBATCH --error=/courses/BINF6430.202530/students/james.ri/logs/salmon_quant_%j.err

shopt -s expand_aliases
CRN=$( basename $( find /courses -maxdepth 1 -type d -name "BINF6430.*" ) | cut -d "." -f2 )

module load singularity/3.10.3
alias salmon="singularity run -B '/courses:/courses,/scratch:/scratch' /courses/BINF6430.${CRN}/shared/singularity_containers/salmon-latest.sif salmon"

BASE=/courses/BINF6430.${CRN}
DATA=${BASE}/data/ReaganData
BAM=${DATA}/aligned_gencode/merged
TRANSOME=${BASE}/shared/gencode/gencode.v46.transcripts.fa

ACC_FILE=${DATA}/merged_replicates.txt
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ACC_FILE})
SAMPLE_SUFFIX="_merged.aligned.toTranscriptome.out.bam"

OUT_DIR="/courses/BINF6430.202530/students/james.ri/STAR_salmon"
mkdir -p ${OUT_DIR}

salmon quant --threads ${SLURM_CPUS_PER_TASK} \
    --targets ${TRANSOME} \
    --numBootstraps 50 \
    --gencode \
    --libType A \
    --seqBias --gcBias --posBias \
    --output ${OUT_DIR}/${SAMPLE} \
    --alignments ${BAM}/${SAMPLE}${SAMPLE_SUFFIX}

