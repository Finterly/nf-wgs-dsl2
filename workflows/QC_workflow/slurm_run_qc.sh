#!/bin/bash
# --- Start of slurm commands -----------
# Request an hour of runtime:
#SBATCH --time=24:00:00
# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH --mem=10G
# Specify a job name:
#SBATCH -J qc
# Specify an output file
# %j is a special variable that is replaced by the JobID when 
# job starts
#SBATCH -o qc-%j.out
#SBATCH -e qc-%j.out
#----- End of slurm commands ----

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                         # e.g. "did my job exceed its memory request?

conda activate bioinfo

INPUT=/users/fhu18/data/shared/S3TIMULATE_2/test_data
OUTPUT=/users/fhu18/data/shared/S3TIMULATE_2/test_results
TRIM=/users/fhu18/nf-wgs-dsl2/workflows/refs/adapters/NexteraPE-custom.fa

nextflow run qc_workflow.nf -profile slurm,apptainer --inputdir $INPUT --outdir $OUTPUT --trimadapter $TRIM 

exit 0