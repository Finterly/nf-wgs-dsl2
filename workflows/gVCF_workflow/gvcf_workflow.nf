#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// nextflow run gvcsf_workflow.nf -profile docker --inputdir /Users/fihu/Documents/GitHub/s3timulate/nf-wgs-dsl2_gVCFs/data/\*

// VARIANT CAlLING PLASMODIUM FALCIPARUM WGS PROCESSING

params.refdir = "$projectDir/../refs/genomes"

log.info """\
    VARIANT CAlLING PLASMODIUM FALCIPARUM WGS - N F   P I P E L I N E
    ===================================
    refdir		  	: ${params.refdir}
	outdir			: ${params.outdir}
    """
    .stripIndent()

// Variant Calling starts here
// Running HaplotypeCaller to generate gVCFs
process variant_calling {
	
	tag "g variant calling ${pair_id}"

	publishDir "${params.outdir}/$pair_id"
       
    input:
	tuple val(pair_id), path(pf_bam)
    path refdir

    output:
    tuple val(pair_id), path("${pair_id}.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14}.g.vcf"), 
    path("${pair_id}.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14}.g.vcf.idx")  

    script:
    """    
	samtools index -bc ${pf_bam}

	for i in 1 2 3  4 5 6 7 8 9 10 11 12 13 14
	    do
	        gatk --java-options "-Xmx${params.gatk_memory}g" HaplotypeCaller -R $refdir/Pf3D7.fasta -I ${pf_bam} -ERC GVCF -ploidy 2 \
	        --native-pair-hmm-threads 16 -O  ${pair_id}.chr"\$i".g.vcf --assembly-region-padding 100 \
	        --max-num-haplotypes-in-population 128 --kmer-size 10 --kmer-size 25 \
	        --min-dangling-branch-length 4 --heterozygosity 0.0029 --indel-heterozygosity 0.0017 --min-assembly-region-size 100 -L $refdir/core_chr"\$i".list -mbq 5 -DF MappingQualityReadFilter --base-quality-score-threshold 12
	    done
    """	
}

workflow.onComplete { 
	println ( workflow.success ? "\nDone!": "Oops .. something went wrong" )
}


workflow {
	/*
	Create 'input_ch' channel that emits for each bam/read pair a
	tuple containing 2 elements: pair_id, bam_pf
	*/
	Channel
        .fromPath(params.input, checkIfExists: true)
		.map {tuple( it.name.split('.sorted')[0], it )}
		.set{input_ch}
		//.ifEmpty{error "Cannot find any reads matching: ${params.reads}"}

	// variant calling
	var_ch = variant_calling(input_ch, params.refdir) 

}
