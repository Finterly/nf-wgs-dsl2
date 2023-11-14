# Plasmodium Falciparum WGS Pipeline (Nextflow DSL 2)

Adapted from: 
- https://github.com/Karaniare/Optimized_GATK4_pipeline (shell script)
- https://github.com/jhoneycuttr/nf-wgs (Nextflow DSL 1)

## Directory Organization
- `data`: suggested folder for raw fastq.gz files
- `results`: suggested folder for output
- `workflows` 
  - `qc.nf`: qc workflow for wgs pipeline
  - `gvcf.nf`: qVCFs workflow for wgs pipeline
- `config`
  - `Apptainer`: the file used to build nf-wgs-dsl2.sif  
  - `Dockerfile`: dockerfile for when running localling 
  - `base.config`: base config file 
  - `envs`: conda envs (under construction :construction:)
- `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
  - `adapters`: folder containing trimmomatic adapter files
  - `genomes`: reference genome files and more
  - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`

# About nextflow.config
|Parameters|Description|
|---|---|
|qc_only|If enabled, only QC workflow is run (default 'false')|
|gvcf_only|If enabled, only gVCF workflow is run (default 'false')|
|inputdir|The folder that contains the input files (default 'data')|
|outputdir|The folder where you want the resulting data to be save (default 'results/results')|
|trimadapter|The adapter used for initial trimming of reads (default 'NexteraPE-custom.fa')|

Additionally, the nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on. The different profiles are listed below, including any setup that is required.

# Parameters in main.nf
|Other Parameters|Description|
|---|---|
|refdir|Reference directory is assumed to be located here '$projectDir/refs/genomes'|
|rscript|Rscript for run quality report is assumed to be located here "$projectDir/refs/run_quality_report.Rmd'|
|reads|If running the pipeline starting from QC, inputdir is assumed to contain the raw reads fastq.gz files ('${params.inputdir}/*_R{1,2}*.fastq.gz')|
|bams|If running the pipeline starting from gVCF, inputdir is assumed to contain the pf bam files ('${params.inputdir}/*/*.sorted.dup.pf.{bam,bam.csi}')|

## Running the Workflow

There are several options for running the Workflow. 

### 1A. Apptainer
Run  `main.nf` on Wynton using Apptainer(singularity container). 

If the apptainer image is not already available, please run the command below to generate the apptainer image. Use `sudo` if necessary.
```bash
apptainer build nf-wgs-dsl2.sif Apptainer
```

And then include the `apptainer` profile on the command line. *Note: you should also include executor you wish to run*

```bash
nextflow run main.nf -profile sge,apptainer
```

Below is an example using some parameter. Please be sure to specify **full paths**!:

```bash
nextflow run main.nf -profile sge,apptainer --inputdir path/input_directory --outdir path/output_directory --trimadapter path/adapters/NexteraPE-custom.fa --qc_only
```

### 1B. Apptainer + Job script
Submit `run_wgs.sh` script as Wynton job. This option is essentially the same as Option 1,  but packaged into a script. 

The `run_wgs.sh` script contains a bash command for running the nextflow workflow using Apptainer. 
You must specify the **full path** to the desired input directory, output directory, and trimmomatic adapter (optional)

snippet from `run_wgs.sh`: 
```
...

INPUT=/path_to/WGS_pipeline_nextflow/data
OUTPUT=/path_to/WGS_pipeline_nextflow/results
TRIM=/wynton/scratch/finterly_WGS_pipeline/workflows/refs/adapters/NexteraPE-custom.fa

nextflow run qc_workflow.nf -profile sge,apptainer --inputdir $INPUT --outdir $OUTPUT --trimadapter $TRIM 
```

`run_qc_workflow.sh` can be run using the following command on Wynton from a log or dev node: 
```bash
qsub -cwd run_qc_workflow.sh
```
You can monitor job progress using `qstat` or viewing the log `cat run_qc_workflow.sh.o#######`. 


### 1C. Docker
Run  `qc_workflow.nf` on local computer using Docker image. The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t finterly/nf-wgs-dsl2 .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
NXF_VER=22.11.0-edge nextflow run qc_workflow.nf -profile docker
```


### Additional information 

- Running nextflow with conda profile under construction :construction: (We prefer Apptainer and Docker anyways...)

- If you need to resume after some processes were successfully executed, add -resume at the end of it
