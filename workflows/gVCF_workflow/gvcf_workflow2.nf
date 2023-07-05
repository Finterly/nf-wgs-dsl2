#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// GENOMIC VARIANT CAlLING PLASMODIUM FALCIPARUM WGS PIPELINE

params.refdir = "$projectDir/../refs/genomes"
params.chrom_range = '1..3'

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
	
	tag "g variant calling ${pair_id}"

	publishDir "${params.outdir}/$chrom"
       
    input:
	tuple val(pair_id), path(pf_bam), path(pf_bam_index)
    path refdir
	val chrom
	

    output:
    tuple val(pair_id), path("${pair_id}.chr${chrom}.g.vcf"), 
    path("${pair_id}.chr${chrom}.g.vcf.idx")  

    script:
    """    
	gatk --java-options "-Xmx${params.gatk_memory}g" HaplotypeCaller -R $refdir/Pf3D7.fasta -I ${pf_bam} -ERC GVCF -ploidy 2 \
	--native-pair-hmm-threads 16 \
	-O  ${pair_id}.chr${chrom}.g.vcf \
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
	--base-quality-score-threshold 12
    """	
}


workflow.onComplete { 
	println ( workflow.success ? "\nDone!": "Oops .. something went wrong" )
}


workflow {
	/*
	Create 'input_ch' channel that emits for each bam/read pair a
	tuple containing 2 elements: pair_id, pf_bam, pf_bam_index
	*/

	Channel
		.fromPath(params.input, checkIfExists: true)
        .filter{it.name.endsWith('.bam')}
        .map{
            tuple(it.name.split('.sorted')[0], it, it + '.csi')
        }
        .set{input_ch}
		//.ifEmpty{error "Cannot find any reads matching: ${params.reads}"}

	input_ch.view()
	
    // Loop over for chromosomes 1 through 14 (default) 
    chroms = set(params.chrom_range.split(','))

	// variant calling
    for (chrom in chroms) {
        g_variant_calling(input_ch, params.refdir, chrom)
    } 
}
