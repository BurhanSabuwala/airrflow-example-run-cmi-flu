#!/bin/bash
#SBATCH --job-name=cmi-flu-tcr
#SBATCH --out="slurm-%j.out"
#SBATCH --error="slurm-%j.err"
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# module load Java/17.0.4
export SQUEUE_FORMAT="%.12i %.8P %.60j %.6u %.2t %.12M %.12l %.24R %.4D %.4C %.6m %.8b %.8f %.10Q"
export NXF_WRAPPER_STAGE_FILE_THRESHOLD='40000'
#export TOWER_WORKSPACE_ID=118983227025632

nextflow pull nf-core/airrflow -r 5.1.0

nextflow run nf-core/airrflow -r 5.1.0 \
-profile singularity \
--mode assembled \
-c immcantation.config \
--input airrflow_tcr_samplesheet_kleinstein.tsv \
--outdir output3 \
--clonal_threshold 0 \
-w work \
-resume
