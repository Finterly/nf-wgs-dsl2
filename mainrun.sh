#!/usr/bin/env bash

REPO=https://github.com/Finterly/nf-wgs-dsl2.git
BRANCH=changes-for-v0.0.6
STAMP=$(date +%H%m%M%S)
MAIN=$BRANCH-$STAMP
git clone -b $BRANCH --single-branch $REPO $MAIN

echo $MAIN
cd $MAIN

SUBMIT="script.sh"
#--------Begin SGE Script-----------#
cat > $SUBMIT <<- PROPERTIES
#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=1G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=24:00:00   # job requires up to 24 hours of runtime
#$ -r n               # if job crashes, it should be restarted

date
hostname

## End-of-job summary, if running as a job
[[ -n "\$JOB_ID" ]] && qstat -j "\$JOB_ID"  # This is useful for debugging and usage purposes,
                                         # e.g. "did my job exceed its memory request?

TRIM=/wynton/home/eppicenter/finterly/kn_test/workflow/adapters/NexteraPE-custom.fa

TIME=time.txt

NXF_VER=22.11.0-edge /usr/bin/time -v -o \$TIME nextflow -q run main.nf --trimadapter \$TRIM -profile sge,apptainer

exit 0

PROPERTIES
#----------End SGE Script----------#

qsub -cwd $SUBMIT

cd ..

exit 0