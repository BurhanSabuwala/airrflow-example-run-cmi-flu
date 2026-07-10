#!/bin/bash
#SBATCH --job-name=cmi-reference-tcr
#SBATCH --out="logs/subj5-slurm-%j.out"
#SBATCH --err="logs/subj5-slurm-%j.err"
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=98G
#SBATCH --partition=day

module load Java/17.0.4
export SQUEUE_FORMAT="%.12i %.8P %.60j %.6u %.2t %.12M %.12l %.24R %.4D %.4C %.6m %.8b %.8f %.10Q"
export NXF_WRAPPER_STAGE_FILE_THRESHOLD='40000'
# export TOWER_WORKSPACE_ID=118983227025632

# nextflow pull nf-core/airrflow -r 5.0.0 -latest

nextflow run nf-core/airrflow -r 5.0.0 \
-profile singularity \
--mode assembled \
--input subject-specific-script/airrflow_bulk_samplesheet_7742.tsv \
--outdir output_per_subject/s7742 \
--reference_fasta References/7742-LI-1/custom_imgtdb.zip \
--clonal_threshold 0 \
-w work-ref-s7742 
