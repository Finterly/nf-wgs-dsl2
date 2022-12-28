#!/usr/bin/env nextflow
 nextflow.enable.dsl=2

// kniare part I of WGS processing for variant calling

ref = file(params.ref)
refdir = file(params.refdir)
rscript = file(params.rscript)
trimadapter = file(params.trimadapter)


log.info """\
    GATK4 OPTIMIZED PART I WGS - N F   P I P E L I N E
    ===================================
    ref           	: ${params.ref}
    refdir		  	: ${params.refdir}
	rscript		  	: ${params.rscript}
    trimadapter		: ${params.trimadapter}
    """
    .stripIndent()


//trimmomatic read trimming
process trimreads {
	
	tag "trim ${pair_id}"	

	publishDir "${params.outdir}/$pair_id",
		saveAs: {filename ->
			if (filename.indexOf("_paired.fq.gz") > 0) "trimmed_pairs/$filename"
			else if (filename.indexOf("_unpaired.fq.gz") > 0) "unpaired/$filename"
			else filename
	}
		
	input:
	tuple val(pair_id), path(reads)

	output:
	tuple val(pair_id), path("trimmed_${pair_id}_R{1,2}_paired.fq.gz"),
	path("trimmed_${pair_id}_R{1,2}_unpaired.fq.gz")

	script:
	"""
	trimmomatic PE ${reads[0]} ${reads[1]} \
	"trimmed_${pair_id}_R1_paired.fq.gz" "trimmed_${pair_id}_R1_unpaired.fq.gz" \
	"trimmed_${pair_id}_R2_paired.fq.gz" "trimmed_${pair_id}_R2_unpaired.fq.gz" \
	ILLUMINACLIP:$trimadapter:2:30:10 LEADING:3 TRAILING:3 MINLEN:3 SLIDINGWINDOW:5:20 -threads ${params.max_threads}
	"""
}

//fastqc on each trimmed read pair
process fastqc {
    
    tag "FASTQC on ${pair_id}"

    publishDir "${params.outdir}/${pair_id}/fastqc"

    input:
    tuple val(pair_id), path(paired_reads), path(unpaired_reads)

    output:
    tuple val(pair_id), path("fastqc_${pair_id}")

    conda 'bioconda::fastqc'
   
    script:
    """
    mkdir fastqc_${pair_id}
    fastqc -o fastqc_${pair_id} -q ${paired_reads}
    """  
}  

//multiqc report
process multiqc {
	
	tag "multiqc on all trimmed_fastqs"

    publishDir params.outdir, mode:'copy'
       
    input:
    tuple val(pair_id), path(fastqc_results)
    
    output:
    file('multiqc_report.html')  
    
    conda 'bioconda::multiqc'

    script:
    """    
    multiqc .
    """
}


// bwa alignment
process bwa_align {
	
	tag "align ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(paired_reads), path(unpaired_reads)

	output:
	tuple val(pair_id), path("${pair_id}.sam")
	
	script:

	"""
	# load modules
	# module load CBI bwa

	# alignment and populate read group header
	bwa mem -t ${params.max_threads} -M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:illumina\\tSM:${pair_id}\\tPU:${pair_id}" $ref ${paired_reads} > ${pair_id}.sam
	"""
}

// sam file sorting
process sam_sort {
	
	tag "sam sorting ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(sam)

	output:
	tuple val(pair_id), path("${pair_id}.sorted.dup.bam"),
	path("${pair_id}.sorted.bam"),
	path("${pair_id}.bam"),
	path("${pair_id}.clean.bam"),
	path("${pair_id}_dup_metrics.txt")
	//file("${pair_id}.sorted.bam.bai")
	//file("${pair_id}.sam")

	"""
	# sam file sorting
	gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" SamFormatConverter -R $ref -I ${pair_id}.sam -O ${pair_id}.bam
    gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" CleanSam -R $ref -I ${pair_id}.bam -O ${pair_id}.clean.bam
    gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" SortSam -R $ref -I ${pair_id}.clean.bam -O ${pair_id}.sorted.bam -SO coordinate --CREATE_INDEX true --TMP_DIR ${PWD}/tmp
    gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" MarkDuplicates -R $ref -I ${pair_id}.sorted.bam -O ${pair_id}.sorted.dup.bam -M ${pair_id}_dup_metrics.txt -ASO coordinate --TMP_DIR ${PWD}/tmp
    """
}

// samtools sorting Pf and human reads
process sort_pf_human {
	
	tag "sort PfHs ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(pf_bam)

	output:
	tuple val(pair_id), 
	path("${pair_id}.sorted.dup.pf.bam"),
	path("${pair_id}.sorted.dup.hs.bam")
	//file("${pair_id}.sorted.dup.pf.bam.csi")

	script:
	"""
	# sorting of Pf and Hs aligned reads
	samtools view -b -h ${pair_id}.sorted.dup.bam -T $ref -L $refdir/Pf3D7_core.bed > ${pair_id}.sorted.dup.pf.bam
	samtools view -b -h ${pair_id}.sorted.dup.bam -T $ref -L $refdir/human.bed > ${pair_id}.sorted.dup.hs.bam
	
	# samtools index -bc ${pair_id}.sorted.dup.pf.bam
	
	# rm ${pair_id}.sorted.dup.bam
	"""	
}


// distribution of Pf read depth by chromosome
process pf_read_depth {
	
	tag "read depth Pf chroms ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir/by_chrom"

	input:
	tuple val(pair_id), path(pf_bam)

	output:
	file("ReadCoverage_final_${pair_id}.tsv")

	script:
	"""
	# load modules
	# module load CBI samtools gatk/4.2.2.0

	samtools index -bc ${pf_bam}

	for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
	    do
	       gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" DepthOfCoverage \
		   -R "$refdir/Pf3D7.fasta" \
		   -O chr"\$i" \
		   -L Pf3D7_"\$i"_v3 \
		   --omit-locus-table true \
		   -I ${pf_bam} --tmp-dir ${PWD}/tmp
	       awk -F"," -v OFS="\t" '{ print \$0, \$(NF+1) = '"chr\$i"' }' chr"\$i".sample_summary > chr"\$i".sample2_summary
	    done

	cat *.sample2_summary | awk '!/sample_id/ {print \$0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15, chromosome' > ReadCoverage_final_${pair_id}.tsv

	"""
}

/*
// pf read depth summary
process pf_read_depth_summary {
	
	tag "read depth summary ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"

	input:
	set pair_id, file(chrom_summary) from read_cover_ch.collect()

	output:
	file("ReadCoverage_final_${pair_id}.tsv")

	script:
	"""
	#!/usr/bin/env bash

	cat $chrom_summary | awk '!/sample_id/ {print \$0}' | sed '1isample_id, total, mean, third_quartile, median, first_quartile, bases_perc_above_15, chromosome' > ReadCoverage_final_${pair_id}.tsv
	"""
}
*/

// insert size calculation
process insert_sizes {
	
	tag "insert sizes ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"

	input:
	tuple val(pair_id), path(pf_bam)

	output:
	tuple val(pair_id), path("${pair_id}.insert.txt"), path("${pair_id}.insert2.txt")
	
	script:
	"""
	# load modules
	# module load CBI gatk/4.2.2.0

	gatk CollectInsertSizeMetrics -I ${pf_bam} -O ${pair_id}.insert.txt -H ${pair_id}_histo.pdf -M 0.05
	awk 'FNR>=8 && FNR<=8 {print \$1,\$3,\$4,\$5,\$6,\$7,\$8,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$NF="${pair_id}"}' ${pair_id}.insert.txt > ${pair_id}.insert2.txt
	
	#rm ${pair_id}.insert.txt 
	"""
}

// insert summary
process insert_summary {

	tag "insert summary"

	publishDir params.outdir, mode:'copy'

	input:
	tuple val(pair_id), path(insert2_collection)

	output:
	tuple val(pair_id), path('InsertSize_Final.tsv')

	script:
	"""
	cat $insert2_collection > InsertSize_Final.tsv
	"""
}

// Pf bam statistics by sample
process pf_bam_stat_per_sample {
	
	tag "Pf bam stat ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"
	
	input: 
	tuple val(pair_id), path(pf_bam)

	output:
	tuple val(pair_id), path("${pair_id}_bamstat_pf_final.tsv"), path("${pair_id}_bamstat_pf.tsv")

	conda 'bioconda::samtools bioconda::datamash'

	script:
	"""
    samtools stats ${pf_bam} | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_pf.tsv
    
    datamash transpose < ${pair_id}_bamstat_pf.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_pf_final.tsv

    #rm ${pair_id}_bamstat_pf.tsv
    """
}

// Pf bam statistic summary
process pf_stat_summary {
	
	tag "Pf stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	tuple val(pair_id), path(bamstat_pf) 

	output:
	tuple val(pair_id), path('Bam_stats_pf_Final.tsv') 

	script:
	"""
	cat $bamstat_pf  | sed '1irow_total_reads_pf	filtered_reads_pf	sequences_pf	is_sorted_pf	1st_fragments_pf	last_fragments_pf	reads_mapped_pf	reads_mapped_and_paired_pf	reads_unmapped_pf	reads_properly_paired_pf	reads_paired_pf	reads_duplicated_pf	reads_MQ0_pf	reads_QC_failed_pf	non_primary_alignments_pf	total_length_pf	total_first_fragment_length_pf	total_last_fragment_length_pf	bases_mapped_pf	bases_mapped_(cigar)_pf	bases_trimmed_pf	bases_duplicated_pf	mismatches_pf	error_rate_pf	average_length_pf	average_first_fragment_length_pf	average_last_fragment_length_pf	maximum_length_pf	maximum_first_fragment_length_pf	maximum_last_fragment_length_pf	average_quality_pf	insert_size_average_pf	insert_size_standard_deviation_pf	inward_oriented pairs_pf	outward_oriented_pairs_pf	pairs_with_other_orientation_pf	pairs_on_different_chromosomes_pf	percentage_of_properly_paired_reads_(%)_pf	sample_id' > Bam_stats_pf_Final.tsv
	#rm *_bamstat_pf_final.tsv
	"""
}

// Hs bam statistics by sample
process hs_bam_stat_per_sample {

	tag "Hs bam stat ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"
	
	input:
	tuple val(pair_id), path(hs_bam)

	output:
	tuple val(pair_id), path("${pair_id}_bamstat_hs_final.tsv"), path("${pair_id}_bamstat_hs.tsv")

	conda 'bioconda::samtools bioconda::datamash'

	script:
	"""
    samtools stats ${hs_bam} | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_hs.tsv

    datamash transpose < ${pair_id}_bamstat_hs.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_hs_final.tsv

    #rm ${pair_id}_bamstat_hs.tsv
    """
}


// Hs bam statistic summary
process hs_stat_summary {
	tag "Hs stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	tuple val(pair_id), path(bamstat_hs)

	output:
	tuple val(pair_id), path('Bam_stats_hs_Final.tsv')

	script:
	"""
	cat $bamstat_hs | sed '1irow_total_reads_hs	filtered_reads_hs	sequences_hs	is_sorted_hs	1st_fragments_hs	last_fragments_hs	reads_mapped_hs	reads_mapped_and_paired_hs	reads_unmapped_hs	reads_properly_paired_hs	reads_paired_hs	reads_duplicated_hs	reads_MQ0_hs	reads_QC_failed_hs	non_primary_alignments_hs	total_length_hs	total_first_fragment_length_hs	total_last_fragment_length_hs	bases_mapped_hs	bases_mapped_(cigar)_hs	bases_trimmed_hs	bases_duplicated_hs	mismatches_hs	error_rate_hs	average_length_hs	average_first_fragment_length_hs	average_last_fragment_length_hs	maximum_length_hs	maximum_first_fragment_length_hs	maximum_last_fragment_length_hs	average_quality_hs	insert_size_average_hs	insert_size_standard_deviation_hs	inward_oriented pairs_hs	outward_oriented_pairs_hs	pairs_with_other_orientation_hs	pairs_on_different_chromosomes_hs	percentage_of_properly_paired_reads_(%)_hs	sample_id' > Bam_stats_hs_Final.tsv
	#rm *_bamstat_hs_final.tsv
	"""
}


// Pf:Hs read ratio calculation
process pf_hs_ratio_calc {
	
	tag "Pf:Hs ratio"

	publishDir params.outdir, mode:'copy'
	
	input: 
	tuple val(pair_id), path(bamsum_pf), path(bamsum_hs)

	output:
	file('Ratios_hs_pf_reads.tsv')

	shell:
	'''
	paste !{bamsum_pf} !{bamsum_hs} | awk -v OFS="\t" '!/per/ {print$39,$7,$46}' | awk '{if($2==0) $2=1}1' |  awk '{if($3==0) $3=1}1' | awk -v OFS="\t" '{print$1,$2,$3,($3/$2)}' | sed '1ireads_mapped_pf, reads_mapped_hs, ratio_hs_pf' > Ratios_hs_pf_reads.tsv
	'''

}


// Rmd run quality report generation
process run_report {
 
	tag "Run quality report"

	publishDir params.outdir, mode:'copy'

	input:
	tuple val(pair_id), path('Bam_stats_pf_Final.tsv'), path('Bam_stats_hs_Final.tsv')
	
	output:
	file('run_quality_report.html')

	script:
	"""
	Rscript -e 'rmarkdown::render(input = "$rscript", output_dir = getwd(), params = list(directory = "${params.outdir}"))'
	"""
}


workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the these reports in your browser --> \nmultiqc = $params.outdir/multiqc_report.html\nrun quality report = $params.outdir/reportrmd.html\nnextflow summary = $params.outdir/report.html": "Oops .. something went wrong" )
}


workflow {

	/*
	Create 'read_pairs' channel that emits for each read pair a
	tuple containing 3 elements: pair_id, R1, R2
	*/

	Channel
        .fromFilePairs(params.reads, checkIfExists: true)
		.set{read_pairs_ch}
		//.ifEmpty{error "Cannot find any reads matching: ${params.reads}"}

    // trim reads
    trimmed_reads_ch = trimreads(read_pairs_ch)
    
    // fastqc report 
    fastqc_ch = fastqc(trimmed_reads_ch)
    // multiqc report --
	multiqc(fastqc_ch.collect()) 

    // bwa alignment
    sam_ch = bwa_align(trimmed_reads_ch)

	// sam file sorting
	sam_sort_ch = sam_sort(sam_ch)
	dup_sort_sam_ch = sam_sort_ch.map{T->[T[0],T[1]]}

	// samtools sorting Pf and human reads
	bam_ch = sort_pf_human(dup_sort_sam_ch)
	pf_bam_ch = bam_ch.map{T->[T[0],T[1]]} // select pf bam
	hs_bam_ch = bam_ch.map{T->[T[0],T[2]]} // select hs bam

	// distribution of Pf read depth by chromosome -- 
	pf_read_depth(pf_bam_ch) 

	// insert size calculation
	inserts_ch = insert_sizes(pf_bam_ch) 
	insert2_ch = inserts_ch.map{T->[T[0],T[2]]} // select *.insert2.txt

	// insert summary -- 
	insert_summary(insert2_ch.collect()) 

	// Pf bam statistics by sample
	pf_bamstat_ch = pf_bam_stat_per_sample(pf_bam_ch)
	pf_final_bamstat_ch = pf_bamstat_ch.map{T->[T[0],T[1]]} // select *_bamstat_pf_final.tsv
	// Pf bam statistic summary
	pf_summary_ch = pf_stat_summary(pf_final_bamstat_ch.collect()) 


	// Hs bam statistics by sample
	hs_bamstat_ch = hs_bam_stat_per_sample(hs_bam_ch)
	hs_final_bamstat_ch = hs_bamstat_ch.map{T->[T[0],T[1]]} // select *_bamstat_hs_final.tsv
	// Hs bam statistic summary
	hs_summary_ch = hs_stat_summary(hs_final_bamstat_ch.collect())

	// Pf:Hs read ratio calculation -- 
	summary_ch = pf_summary_ch.join(hs_summary_ch) //join 
	pf_hs_ratio_calc(summary_ch) 

	// Rmd run quality report generation -- 
	pf_summary_collect_ch = pf_summary_ch.collect() //collect
	hs_summary_collect_ch = hs_summary_ch.collect() //collect
	summary_collect_ch = pf_summary_collect_ch.join(hs_summary_collect_ch) //join 
	run_report(summary_collect_ch)

}