#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// GENOMIC VARIANT CAlLING PLASMODIUM FALCIPARUM WGS PIPELINE

params.refdir = "$projectDir/../refs/genomes"

log.info """\
    GENOMIC VARIANT CAlLING PLASMODIUM FALCIPARUM WGS - N F   P I P E L I N E
    ===================================
    refdir		  	: ${params.refdir}
	outdir			: ${params.outdir}
    """
    .stripIndent()

// Genomic Variant Calling starts here
// Running HaplotypeCaller to generate gVCFs
process g_variant_calling {
	
	tag "g variant calling ${pair_id} chrom ${chrom}"

    publishDir "${params.outdir}/${chrom}", mode:'copy'
       
    input:
	tuple val(pair_id), path(pf_bam_and_idx), val(chrom)
    path refdir

    output:
    tuple path("${pair_id}.chr${chrom}.g.vcf"), path("${pair_id}.chr${chrom}.g.vcf.idx")  

    script:
    """    
	gatk --java-options "-Xmx${params.gatk_memory}g" HaplotypeCaller \
	-R $refdir/Pf3D7.fasta \
	-I ${pf_bam_and_idx[0]} \
	-ERC GVCF \
	-ploidy 2 \
	--native-pair-hmm-threads 16 \
	-O ${pair_id}.chr${chrom}.g.vcf \
	--assembly-region-padding 100 \
	--max-num-haplotypes-in-population 128 \
	--kmer-size 10 \
	--kmer-size 25 \
	--min-dangling-branch-length 4 \
	--heterozygosity 0.0029 \
	--indel-heterozygosity 0.0017 \
	--min-assembly-region-size 100 \
	-L $refdir/core_chr${chrom}.list \
	-mbq 5 \
	-DF MappingQualityReadFilter \
	--base-quality-score-threshold 122
    """
}


workflow.onComplete { 
	println ( workflow.success ? "\nDone!": "Oops .. something went wrong" )
}

workflow {

	// Create 'bam_ch' channel that emits for each sample a tuple containing 3 elements: pair_id, pf_bam, pf_bam_index
	bam_ch = Channel.fromFilePairs(params.input, checkIfExists: true)
	// Create 'chrom_ch' channel which emits a number for pf chromosomes 1 through 14
	chrom_ch = Channel.from(1,2,3,4,5,6)

	// Combine 'bam_ch' and 'chrom_ch' to create 'input_ch' channel that emits for each sample
	// a tuple containing 4 elements: pair_id, pf_bam, pf_bam_index, chromosome number
	input_ch = bam_ch.combine(chrom_ch)
	
	// variant calling
	var_ch = g_variant_calling(input_ch, params.refdir) 
	
}
