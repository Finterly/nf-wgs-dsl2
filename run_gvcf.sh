#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=16G     # job requires up to 16 GiB of RAM per slot
#$ -l scratch=60G      # job requires up to 20 GiB of local /scratch space
#$ -l h_rt=120:00:00   # job requires up to 24 hours of runtime
#$ -r n               # if job crashes, it should be restarted

date
hostname

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  
# This is useful for debugging and usage purposes,
# e.g. "did my job exceed its memory request?

module load CBI miniconda3/23.5.2-0-py311 #wynton 
conda activate nextflow_env #optional

# Make sure to modify Paths! 
WORKDIR=/wynton/scratch/finterly_temp/run_FABG01_gvcf/work_gvcf
INPUT=/wynton/scratch/finterly_temp/run_FABG01_qc/results/final_bams
OUTPUT=/wynton/scratch/finterly_temp/run_FABG01_gvcf/results_gvcf

# Nextflow pipeline stored in shared drive
nextflow run /wynton/home/eppicenter/shared/WGS_pipeline_nextflow/main.nf \
-profile sge,apptainer \
-work-dir $WORKDIR \
--inputdir $INPUT \
--outputdir $OUTPUT \
--qc_only true

exit 0