#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=16G     # job requires up to 16 GiB of RAM per slot
#$ -l scratch=20G      # job requires up to 20 GiB of local /scratch space
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime
#$ -r n               # if job crashes, it should be restarted

date
hostname

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  
# This is useful for debugging and usage purposes,
# e.g. "did my job exceed its memory request?

conda activate nextflow_env #optional

# Make sure to modify Paths! 
INPUT=/path_to/WGS_pipeline_nextflow/data
OUTPUT=/path_to/WGS_pipeline_nextflow/results

nextflow run main.nf -profile sge,apptainer \
--inputdir $INPUT \
--outputdir $OUTPUT

exit 0