#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Genomic Variant Calling Workflow

// Running HaplotypeCaller to generate gVCFs
process g_variant_calling {
    
    tag "g variant calling ${pair_id} chr${chrom}"
    label 'big_mem'
    
    publishDir "${params.outputdir}/chr${chrom}", mode:'copy'
       
    input:
    tuple val(pair_id), path(pf_bam), path(pf_bam_index), val(chrom)
    path genomes_dir

    output:
    tuple path("${pair_id}.chr${chrom}.g.vcf"), path("${pair_id}.chr${chrom}.g.vcf.idx"), path("${pair_id}.chr${chrom}_log.txt")  

    script:
    """    
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \
    -R $genomes_dir/Pf3D7.fasta \
    -I ${pf_bam} \
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
    -L $genomes_dir/core_chr${chrom}.list \
    -mbq 5 \
    -DF MappingQualityReadFilter \
    --base-quality-score-threshold 12 \
    > ${pair_id}.chr${chrom}_log.txt 2>&1
    """
}


workflow.onComplete { 
    println ( workflow.success ? "\ngVCF run complete!": "Oops .. something went wrong" )
}

workflow GVCF {

    take: pf_bam_ch
    main: 
        // Create 'chrom_ch' channel which emits a number for pf chromosomes 1 through 14
        chrom_ch = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14)

        // Combine 'bam_ch' and 'chrom_ch' to create 'input_ch' channel that emits for each sample
        // a tuple containing 4 elements: pair_id, pf_bam, pf_bam_index, chromosome number
        input_ch = pf_bam_ch.combine(chrom_ch)
    
        // variant calling
        var_ch = g_variant_calling(input_ch, params.genomes_dir) 

    emit: var_ch
    
}
