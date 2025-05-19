#!/bin/bash
#SBATCH --job-name=FastQC_Run                                  # Job name
#SBATCH --partition=courses                                   # Use course partition
#SBATCH -N 1                                                 # Number of nodes
#SBATCH -c 8                                                 # Number of CPU threads
#SBATCH --mem=32G                                            # Memory allocation
#SBATCH -t 4:00:00                                           # Time limit (HH:MM:SS)
#SBATCH --mail-type=END,FAIL                                 # Email notifications
#SBATCH --mail-user=james.ri@northeastern.edu            # Replace with your email
#SBATCH --out=/courses/BINF6430.202530/students/${USER}/logs/%x_%j.log   # Log file
#SBATCH --error=/courses/BINF6430.202530/students/${USER}/logs/%x_%j.err # Error log

### VARIABLES ###
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

# Load FastQC module
module load fastqc

# Set directory paths
FASTQ_DIR=/courses/BINF6430.202530/data/ReaganData/Reagan_PE85_TakaraPicoV2_HC_CM_10042022/

# Create output directory
mkdir -p /scratch/${USER}/results/fastqc

# Run FastQC on all FASTQ files in the dataset
fastqc -t 8 -o /scratch/${USER}/results/fastqc ${FASTQ_DIR}/*/*.fastq.gz
