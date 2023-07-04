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
// Combining gvcfs
process combine_gvcfs {
	
	tag "combine gvcfs ${chrom}"

	publishDir "${params.outdir}"
       
    input:
	tuple val(chrom), path(gvcfs)

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
	#    gatk --java-options "-Xmx${params.gatk_memory}g" GenomicsDBImport \
	#    --sample-name-map gvcf_chr"$i"_list.tsv \
	#    --genomicsdb-workspace-path chr"$i"_database \
	#    --intervals $ref_dir/core_chr"$i".list --batch-size 100 \
	#    --reader-threads 24 --genomicsdb-segment-size 8048576 \
	#    --genomicsdb-vcf-buffer-size 160384
	#done
    """   
}


// Running joint genotyping by genomic segment 
// (the core genome of each chromosome is split into smaller genomic regions) via multiple slurm jobs
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
	for i in $chrom
	  do
	   cd $ref_dir
	   for j in $(cat $region)
	     do
	     cd $vcf_dir
	     echo -e "#\!/bin/bash" > genotype_chr"$i"_"$j".sh
	     echo -e "#SBATCH -J Genotype_chr"$i"_"$j"" >> genotype_chr"$i"_"$j".sh
	     echo -e "#SBATCH -t 120:00:00" >> genotype_chr"$i"_"$j".sh
	     echo -e "#SBATCH -c 8" >> genotype_chr"$i"_"$j".sh
	     echo -e "#SBATCH --mem-per-cpu=10g" >> genotype_chr"$i"_"$j".sh
	     echo -e "module load gatk/4.2.2.0 samtools/ bcftools/1.9 bcftools/1.9" >> genotype_chr"$i"_"$j".sh
	     echo -e "ulimit -c unlimited" >> genotype_chr"$i"_"$j".sh
	     echo -e "module load java/jdk-17.0.2" >> genotype_chr"$i"_"$j".sh
	     echo -e "gatk --java-options "-Xmx80g -Xms80g"  GenotypeGVCFs --genomicsdb-use-bcf-codec true  -R $ref_dir/Pf3D7.fasta -V gendb://$vcf_dir/chr"$i"_database --max-genotype-count 1024 -O $vcf_dir/chr"$i"_part"$j".vcf.gz --tmp-dir $ref_dir -stand-call-conf 30 -L "$j"" >> genotype_chr"$i"_"$j".sh
	     chmod +x genotype_chr"$i"_"$j".sh
	     sed -i 's/-Xmx80g -Xms80g/"-Xmx80g -Xms80g"/g' genotype_chr"$i"_"$j".sh
	     sbatch genotype_chr"$i"_"$j".sh
	   done
	done
	}
    """   
}

workflow.onComplete { 
	println ( workflow.success ? "\nDone!": "Oops .. something went wrong" )
}


workflow {
	/*
	Create 'input_ch' channel that emits for each chrom a
	tuple containing multiple elements: chrom
	*/
	Channel
        .fromPath(params.input, checkIfExists: true)
		.map {tuple( it.name.split('.')[1].split('.g.vcf*')[0], it )}
		.set{input_ch}
		//.ifEmpty{error "Cannot find any reads matching: ${params.reads}"}

	// variant calling
	combine_gvcfs(input_ch) 

}
