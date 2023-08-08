#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// VARIANT CAlLING PLASMODIUM FALCIPARUM WGS PIPELINE

params.refdir = "$projectDir/../refs/genomes"

log.info """\
    VARIANT CAlLING PLASMODIUM FALCIPARUM WGS - N F   P I P E L I N E
    ===================================
    refdir		  	: ${params.refdir}
	outdir			: ${params.outdir}
    """
    .stripIndent()


// Combining gVCFs into a genomic database
process combine_gvcfs {
	
	tag "combine gvcfs ${chrom}"

	publishDir "${params.outdir}"
       
    input:
	tuple val(pair_id), path(g_vcf_and_idx)

    output:
    tuple val(chrom)

    script:
    """   
	echo $chrom 
	echo $gvcfs
	#ref_dir="/users/kniare/data/shared/WGS_works_Karamoko/scripts/Data_WGS_UNC"
	#cd $vcf_dir
	#for i in $chrom
	#  do
	
	gatk --java-options "-Xmx${params.gatk_memory}g" GenomicsDBImport \
	--sample-name-map gvcf_chr"$i"_list.tsv \
	--genomicsdb-workspace-path chr"$i"_database \
	--intervals $ref_dir/core_chr"$i".list \
	--batch-size 100 \
	--reader-threads 24 \
	--genomicsdb-segment-size 8048576 \
	--genomicsdb-vcf-buffer-size 160384
	
	#done
    """   
}


// Running joint genotyping by genomic segment (the core genome of each chromosome 
// is split into smaller genomic regions) via multiple slurm jobs
process genotype {
	
	tag "combine gvcfs ${pair_id}"

	publishDir "${params.outdir}/$pair_id"
       
    input:
	tuple val(pair_id), path(pf_bam)
    path refdir

    output:
    tuple val(pair_id), path("${pair_id}.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14}.g.vcf"), 
    path("${pair_id}.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14}.g.vcf.idx")  

    script:
    """   	
	##### 
	#for i in $chrom
	#  do
	#   cd $ref_dir
	#   for j in cat $region
	#     do
	#     cd $vcf_dir
		text_data=$(cat <<EOT
		#!/bin/bash
		#SBATCH -J Genotype_chr"$i"_"$j"
		#SBATCH -t 120:00:00
		#SBATCH -c 8
		#SBATCH --mem-per-cpu=10g
		module load gatk/4.2.2.0 samtools/ bcftools/1.9 bcftools/1.9
		ulimit -c unlimited
		module load java/jdk-17.0.2
		gatk --java-options "-Xmx80g -Xms80g"  GenotypeGVCFs \
		--genomicsdb-use-bcf-codec true  \
		-R $ref_dir/Pf3D7.fasta \
		-V gendb://$vcf_dir/chr"$i"_database \
		--max-genotype-count 1024 \
		-O $vcf_dir/chr"$i"_part"$j".vcf.gz \
		--tmp-dir $ref_dir -stand-call-conf 30 \
		-L "$j"" >> genotype_chr"$i"_"$j".sh 
		EOT
		)
		echo "$text_data" > genotype_chr"$i"_"$j".sh



	#      chmod +x genotype_chr"$i"_"$j".sh
	#     sed -i 's/-Xmx80g -Xms80g/"-Xmx80g -Xms80g"/g' genotype_chr"$i"_"$j".sh
	#     sbatch genotype_chr"$i"_"$j".sh
	#   done
	#done
	#}
    """   


	
}

workflow.onComplete { 
	println ( workflow.success ? "\nDone!": "Oops .. something went wrong" )
}


workflow {

	// Create 'chrom_ch' channel which emits a number for pf chromosomes 1 through 14 (gvcf folders)
	chrom_ch = Channel.from(1,2,3,4,5,6)
	
	// Create 'bam_ch' channel that emits for each sample a tuple containing 3 elements: pair_id, pf_bam, pf_bam_index
	bam_ch = Channel.fromFilePairs(params.input, checkIfExists: true)



	// Create 'input_ch' channel that emits for each sample a tuple containing 3 elements: pair_id, g_vcf, g_vcf_index
	Channel
        .fromFilePairs(params.reads, checkIfExists: true)
		.set{input_ch}

	input_ch.view()
		//.ifEmpty{error "Cannot find any reads matching: ${params.reads}"}

	// variant calling
	combine_gvcfs(input_ch) 

}

