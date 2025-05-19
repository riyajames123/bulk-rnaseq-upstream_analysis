# bulk-rnaseq-upstream_analysis
This repository contains scripts and QC reports used for upstream bulk RNA-seq analysis, including quality control, alignment, quantification, and summary reporting.

## 🚫 Data Note
⚠️ Raw FASTQ files are **not included** due to data confidentiality policies. File paths in scripts are absolute and environment-specific.

## 📂 Structure
- `scripts/`: Shell scripts for each pipeline step
- `reports/`: MultiQC reports from pre- and post-alignment

## 🧪 Tools Used
- **FastQC**: Quality check of raw FASTQ files
- **MultiQC**: Aggregated QC reporting
- **STAR**: Genome alignment and indexing
- **Salmon**: Transcript-level quantification

## 🧬 Pipeline Steps
1. Quality check: `scripts/fastqc_single.sh`
2. Aggregate QC: `multiqc_report.html`
3. STAR genome indexing: `scripts/star_idx_gencode.sh`
4. Alignment: `scripts/star_align.sh`
5. Merging replicates: `scripts/merged_replicates.sh`
6. Quantification with Salmon: `scripts/salmon_quant.sh`
7. Final QC (post-alignment): `rnaseq_aligned_report.html`

## 🔄 Reproducibility
Paths are currently absolute and reflect the original working environment. Modify to suit your setup.
