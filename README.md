# Plasmodium Falciparum WGS Pipeline (Nextflow DSL 2)

Adapted from: 
- https://github.com/Karaniare/Optimized_GATK4_pipeline (shell script)
- https://github.com/jhoneycuttr/nf-wgs (Nextflow DSL 1)

## Directory Organization
- `data`: suggested folder for raw fastq.gz files
- `results`: suggested folder output
- `workflows` 
  - `config`
    - `Apptainer`: the file used to build nf-wgs-dsl2.sif  
    - `Dockerfile`: dockerfile for when running localling 
    - `envs`: conda envs (under construction :construction:)
  - `gVCF_workflow`: qVCFs workflow for wgs pipeline
    - `gvcf_workflow.nf`: nextflow workflow for variant calling  
    - `nextflow.config`: nextflow config file for variant calling workflow
    - `run_gvcf_script.sh`: script for running `gvcf_workflow.nf`
  - `QC_workflow`: qc workflow for wgs pipeline
    - `nextflow.config`: nextflow config file for qc workflow  
    - `qc_workflow.nf`: nextflow workflow for qc   
    - `run_qc_workflow.sh`: script for running `qc_workflow.nf`
  - `refs`: reference files used by both `QC_workflow` and `gVCF_workflow`
    - `adapters`: folder containing trimmomatic adapter files
    - `genomes`: reference genome files and more
    - `run_quality_report.Rmd`: r script for quality report used in `QC_workflow`

# 1. QC Workflow 

About the nextflow.config file:

|Parameters|Description|
|---|---|
|inputdir|The folder that contains all the fastq files (default 'data')|
|outdir|The folder where you want the resulting data to be save (default 'results/qc_results')|
|trimadapter|The adapter used for initial trimming of reads (default 'TruSeq3-PE.fa')|

|Other Parameters|Description|
|---|---|
|reads|The fastq files in the inputdir folder|

Additionally, the nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on. The different profiles are listed below, including any setup that is required.

## Running the QC_workflow

There are several options for running the QC_workflow. 

### 1A. Run  `qc_workflow.nf` on Wynton using Apptainer(singularity container) 
If the apptainer image is not already available, please run the command below to generate the apptainer image. Use `sudo` if necessary.
```bash
apptainer build nf-wgs-dsl2.sif Apptainer
```

And then include the `apptainer` profile on the command line. *Note: you should also include executor you wish to run*

```bash
NXF_VER=22.11.0-edge nextflow run qc_workflow.nf -profile sge,apptainer
```

Below is an example using some parameter. Please be sure to specify **full paths**!:

```bash
NXF_VER=22.11.0-edge nextflow run qc_workflow.nf -profile sge,apptainer --inputdir path/input_directory --outdir path/output_directory --trimadapter path/adapters/NexteraPE-custom.fa

```

### 1B. Submit `run_qc_workflow.sh` script as Wynton job 

This option is essentially the same as Option 1,  but packaged into a script. 

The `run_qc_workflow.sh` script contains a bash command for running the nextflow workflow using Apptainer. 
You must specify the **full path** to the desired input directory, output directory, and trimmomatic adapter (optional)

snippet from `run_qc_workflow.sh`: 
```
...

INPUT=/path_to/WGS_pipeline_nextflow/data
OUTPUT=/path_to/WGS_pipeline_nextflow/results
TRIM=/wynton/scratch/finterly_WGS_pipeline/workflows/refs/adapters/NexteraPE-custom.fa

NXF_VER=22.11.0-edge nextflow run qc_workflow.nf -profile sge,apptainer --inputdir $INPUT --outdir $OUTPUT --trimadapter $TRIM 
```

`run_qc_workflow.sh` can be run using the following command on Wynton from a log or dev node: 
```bash
qsub -cwd run_qc_workflow.sh
```
You can monitor job progress using `qstat` or viewing the log `cat run_qc_workflow.sh.o#######`. 


### 1C. Run  `qc_workflow.nf` on local computer using Docker image 

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t finterly/nf-wgs-dsl2 .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
NXF_VER=22.11.0-edge nextflow run qc_workflow.nf -profile docker
```


# 2. gVCF Workflow 

Notes: We expect to run `gVCF_workflow` following `QC_workflow`, therefore the default input directory for `gVCF_workflow`is the default output directory for `QC_workflow`. 

About the nextflow.config file:

|Parameters|Description|
|---|---|
|inputdir|The folder that contains all the bam files. (default 'results/qc_results')|
|outdir|The folder where you want the resulting data to be save (default 'results/gvcf_results')|

|Other Parameters|Description|
|---|---|
|reads|The bam files in the inputdir folder|

Additionally, the nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on. The different profiles are listed below, including any setup that is required.

## Running the gVCF_workflow

### 2A. Run  `gvcf_workflow.nf` on Wynton using Apptainer(singularity container) 
If the apptainer image is not already available, please run the command below to generate the apptainer image. Use `sudo` if necessary.
```bash
apptainer build nf-wgs-dsl2.sif Apptainer
```

And then include the `apptainer` profile on the command line. *Note: you should also include executor you wish to run*

```bash
NXF_VER=22.11.0-edge nextflow run gvcf_workflow.nf -profile sge,apptainer
```

Below is an example using some parameter. Please be sure to specify **full paths**!:

```bash
NXF_VER=22.11.0-edge nextflow run gvcf_workflow.nf -profile sge,apptainer --inputdir path/input_directory --outdir path/output_directory

```

### 2B. Submit `run_gvcf_workflow.sh` script as Wynton job 

This option is essentially the same as Option 1,  but packaged into a script. 

The `run_gvcf_workflow.sh` script contains a bash command for running the nextflow workflow using Apptainer. 
You must specify the **full path** to the desired input directory, output directory, and trimmomatic adapter (optional)

snippet from `run_gvcf_workflow.sh`: 
```
...

INPUT=/path_to/WGS_pipeline_nextflow/data
OUTPUT=/path_to/WGS_pipeline_nextflow/results

NXF_VER=22.11.0-edge nextflow run gvcf_workflow.nf -profile sge,apptainer --inputdir $INPUT --outdir $OUTPUT
```

`run_gvcf_workflow.sh` can be run using the following command on Wynton from a log or dev node: 
```bash
qsub -cwd run_gvcf_workflow.sh
```
You can monitor job progress using `qstat` or viewing the log `cat run_gvcf_workflow.sh.o#######`. 



### 3C. Run  `gvcf_workflow.nf` on local computer using Docker image 

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t finterly/nf-wgs-dsl2 .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
NXF_VER=22.11.0-edge nextflow run gvcf_workflow.nf -profile docker
```


### Additional information 

- Running nextflow with conda profile under construction :construction: (We prefer Apptainer and Docker anyways...)

- If you need to resume after some processes were successfully executed, add -resume at the end of it
