# Optimized GATK4 Pipeline 
## Part 1 Nextflow DSL 2 Workflow

Adapted from: 
- https://github.com/Karaniare/Optimized_GATK4_pipeline (shell script)
- https://github.com/jhoneycuttr/nf-wgs (Nextflow DSL 1)

## Setup

### Setting Parameters

Modify the nextflow.config file:

|Parameters|Description|
|---|---|
|inputdir|The folder that contains all the fastq files (default 'data')|
|outdir|The folder where you want the resulting data to be save (default 'results')|
|refdir|The folder that contains reference genomes and bed files (default 'genomes')|
|trimadapter|The adapter used for initial trimming of reads (default 'TruSeq3-PE.fa')|

|Other Parameters|Description|
|---|---|
|reads|The fastq files in the inputdir folder|
|ref|The reference genome (default 'Pf3D7_human.fa')|
|rscript|The r script for generating report|

Additionally, the nextflow parameter `-profile` can be use to target the infrastructure you wish to run the pipeline on. The different profiles are listed below, including any setup that is required.
    
### Singularity

Under construction :construction: !

If using singularity, please run the command below to generate the singularity image. Use `sudo` if necessary.

```bash
apptainer build nf-wgs-dsl2.sif Apptainer
```

And then include the `singularity` profile on the command line. 

*Note: you should also include executor you wish to run*

```bash
nextflow run main.nf -profile sge,singularity
```

Below is an example using the genome parameter:

```bash
nextflow run main.nf --readDIR ~/Documents/MAD4HATTER_example_data/single -w ~/Documents/work --target v4 -profile sge,singularity --genome PlasmoDB-59_Pfalciparum3D7_Genome.fasta -c conf/custom.config
```

### Docker

The pipeline can be easily run with docker and is the recommended way to run it when not using an HPC.

Follow the steps below to setup your docker image:

*Note: [docker](https://www.docker.com/) is a prerequisite.*

```bash
docker build -t finterly/nf-wgs-dsl2 .
```

And you're done! To run the pipeline, simply add `-profile docker`. 

```bash
nextflow run main.nf -profile docker
```

### Conda

Under construction :construction: !

To use conda, you must first install either [conda](https://docs.conda.io/en/latest/) or [miniconda](https://docs.conda.io/en/latest/miniconda.html). Once installed, include the `conda` profile on the command line.

```bash
nextflow run main.nf --readDIR single --target v3 -profile conda
```

### Customizing for your institution

There is a file named `custom.config` in `conf/` that can be used to tailor processes to your environment. By default,
this file is used to tailor the pipeline for Wynton HPC at UCSF. This file may be altered to fit your institution's profile.

### Examples 
Under construction :construction: !

Potential ways to execute the pipeline:

```bash
# local executor
nextflow run main.nf

# with sge + singularity and Nextera trimadapter
nextflow run main.nf -profile sge,singularity --trimadapter ./adapters/NexteraPE-PE.fa

# with a profile (currently only supports sge)
nextflow run main.nf -profile sge

```

If you need to resume after some processes were successfully executed, add -resume at the end of it