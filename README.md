# Plasmodium Falciparum WGS Pipeline 
## (Nextflow DSL 2)

Adapted from: 
- https://github.com/Karaniare/Optimized_GATK4_pipeline (shell script)
- https://github.com/jhoneycuttr/nf-wgs (Nextflow DSL 1)

## Overview
- `main.nf`: WGS workflow 
- `nextflow.config`: config file
- `workflows` 
  - `qc.nf`: QC sub-workflow 
  - `gvcf.nf`: GVCF sub-workflow
- `config`
  - `Apptainer`: file used to build nf-wgs-dsl2.sif  
  - `Dockerfile`: file for building docker image 
  - `base.config`: base config file 
  - `envs`: conda envs (under construction :construction:)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`
- *`data`: suggested directory for input files*
- *`results`: suggested directory for output*

## Parameters

### nextflow.config
|Parameters|Description|
|---|---|
|qc_only|If enabled, only QC workflow is run (default 'false')|
|gvcf_only|If enabled, only gVCF workflow is run (default 'false')|
|inputdir|The folder that contains the input files (default 'data')|
|outputdir|The folder where you want the resulting data to be save (default 'results/results')|
|trimadapter|The adapter used for initial trimming of reads (default 'NexteraPE-custom.fa')|

The nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on.

#### Parameters in main.nf
|Parameters|Description|
|---|---|
|refdir|Reference directory is assumed to be located here '$projectDir/refs/genomes'|
|rscript|Rscript for run quality report is assumed to be located here "$projectDir/refs/run_quality_report.Rmd'|
|reads|If running the pipeline starting from QC, inputdir is assumed to contain the raw reads fastq.gz files|
|bams|If running the pipeline starting from gVCF, inputdir is assumed to contain the pf bam files and .csi index files|

## Running the Workflow

There are several options for running the workflow. 

### Apptainer (SGE)
Run  `main.nf` on Wynton using Apptainer(singularity container). 

If the apptainer image is not already available, please run the command below to generate the apptainer image. Use `sudo` if necessary.
```bash
apptainer build nf-wgs-dsl2.sif Apptainer
```

**Run QC then GVCF workflow (default)** 

To run on Wynton, include the `apptainer` profile on the command line and the executor you wish to run (`sge`). 

```bash
nextflow run main.nf -profile sge,apptainer
```

**Run only QC workflow** 

Enable `--qc_only`

Optional: specify full paths for `--inputdir`, `--outputdir`, and `--trimadapter`. 

```bash
nextflow run main.nf \
--qc_only \
-profile sge,apptainer \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory \
--trimadapter path/adapters/NexteraPE-custom.fa 
```

**Run only GVCF workflow** 

Enable `--qc_only`  

Note: for gvcf only runs, `--inputdir` points to the directory containing the `.sorted.dup.pf.bam` and `.sorted.dup.pf.bam.csi` files. 

```bash
nextflow run main.nf \
--gvcf_only \
-profile sge,apptainer \
--inputdir path/input_directory_bam \
--outputdir path/output_directory
```

#### Running a Wynton Job (SGE)
On Wynton, it is recommended to submit runs as jobs. 

Below is an example script `run_wgs.sh` which contains a bash command for running the nextflow workflow using Apptainer. 

```
...
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
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  
# This is useful for debugging and usage purposes,
# e.g. "did my job exceed its memory request?

conda activate nextflow_env #optional

INPUT=/path_to/WGS_pipeline_nextflow/data
OUTPUT=/path_to/WGS_pipeline_nextflow/results

nextflow run main.nf -profile sge,apptainer \
--inputdir $INPUT \
--outputdir $OUTPUT

exit 0
```

`run_main.sh` can be run using the following command on Wynton from a log or dev node: 
```bash
qsub -cwd run_wgs.sh
```
You can monitor job progress using `qstat` or viewing the log `cat run_main.sh.o#######`. 


### Docker (Locally)
Run  `main.nf` on local computer using Docker image. The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t finterly/nf-wgs-dsl2 .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
nextflow run qc_workflow.nf -profile docker
```

**Some Other Examples**

```bash
# QC Only
nextflow run main.nf \
--qc_only \
-profile docker \
--inputdir path/input_directory_fastq \
--outputdir path/output_directory \
--trimadapter path/adapters/NexteraPE-custom.fa 

#GVCF Only
nextflow run main.nf \
--gvcf_only \
-profile docker \
--inputdir path/input_directory_bam \
--outputdir path/output_directory 
```


### Additional information 

- Running nextflow with conda profile under construction :construction: (We prefer Apptainer and Docker anyways...)

- If you need to resume after some processes were successfully executed, add -resume at the end of it
