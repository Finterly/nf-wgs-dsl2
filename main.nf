#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Part I of WGS processing for variant calling

params.refdir = "$projectDir/genomes"
params.rscript = "$projectDir/run_quality_report.Rmd"


log.info """\
    GATK4 OPTIMIZED PART I WGS - N F   P I P E L I N E
    ===================================
    refdir		  	: ${params.refdir}
	rscript		  	: ${params.rscript}
    trimadapter		: ${params.trimadapter}
	outdir			: ${params.outdir}
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
	path trimadapter

	output:
	tuple val(pair_id), path("trimmed_${pair_id}_R{1,2}_paired.fq.gz"),
	path("trimmed_${pair_id}_R{1,2}_unpaired.fq.gz")

	script:
	"""
	trimmomatic PE ${reads[0]} ${reads[1]} \
	"trimmed_${pair_id}_R1_paired.fq.gz" "trimmed_${pair_id}_R1_unpaired.fq.gz" "trimmed_${pair_id}_R2_paired.fq.gz" "trimmed_${pair_id}_R2_unpaired.fq.gz" ILLUMINACLIP:$trimadapter:2:30:10 LEADING:3 TRAILING:3 MINLEN:36 SLIDINGWINDOW:5:20 -threads ${params.max_threads}
	"""
}

//fastqc on each trimmed read pair
process fastqc {
    
    tag "FASTQC on ${pair_id}"

    publishDir "${params.outdir}/${pair_id}/fastqc"

    input:
    tuple val(pair_id), path(paired_reads), path(unpaired_reads)

    output:
    path("fastqc_${pair_id}")

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
    path(fastqc_results)
    
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
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.sam")
	
	script:
	"""
	bwa mem -t ${params.max_threads} \
	-M -R "@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:illumina\\tSM:${pair_id}\\tPU:${pair_id}" \
	$refdir/Pf3D7_human.fa ${paired_reads} > ${pair_id}.sam
	"""
}

// sam format converter
process sam_format_converter {
	
	tag "sam format converter ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(sam_file)
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.bam")

	script:
	"""
	gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" SamFormatConverter \
	-R $refdir/Pf3D7_human.fa \
	-I ${sam_file} \
	-O ${pair_id}.bam
    
	# rm ${sam_file}
	"""
}

// sam clean
process sam_clean {
	
	tag "sam cleaning ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(bam_file)
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.clean.bam")

	script:
	"""
    gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" CleanSam \
	-R $refdir/Pf3D7_human.fa \
	-I ${bam_file} \
	-O ${pair_id}.clean.bam
    
	# rm ${bam_file}
	"""
}

// sam file sorting
process sam_sort {
	
	tag "sam sorting ${pair_id}"
	label 'big_mem'
	scratch true

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(clean_bam)
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.sorted.bam")

	afterScript "rm -rf TMP"

	script: 
	"""
	mkdir -p TMP

	# sam file sorting
    gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" SortSam \
	-R $refdir/Pf3D7_human.fa \
	-I ${clean_bam} \
	-O ${pair_id}.sorted.bam \
	-SO coordinate \
	--CREATE_INDEX true \
	--TMP_DIR TMP

   	# rm ${clean_bam}
    """
}

// mark duplicates
process sam_duplicates {
	
	tag "sam mark duplicates ${pair_id}"
	label 'big_mem'
	scratch true

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(sorted_bam)
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.sorted.dup.bam"),
	path("${pair_id}_dup_metrics.txt")

	afterScript "rm -rf TMP"

	script:
	"""
	mkdir -p TMP

    gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" MarkDuplicates \
	-R $refdir/Pf3D7_human.fa \
	-I ${sorted_bam} \
	-O ${pair_id}.sorted.dup.bam \
	-M ${pair_id}_dup_metrics.txt \
	-ASO coordinate \
	--TMP_DIR TMP
    
	# rm ${sorted_bam}
    """
}
// subsampling Plasmodium falciparum specific bams
process target_pf {
	
	tag "target Pf ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(sorted_dup_bam), path(dup_metrics_txt)
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.sorted.dup.pf.bam")

	script:
	"""
	samtools view -b -h ${sorted_dup_bam} \
	-T $refdir/Pf3D7.fasta \
	-L $refdir/Pf3D7_core.bed > ${pair_id}.sorted.dup.pf.bam
	"""	
}

// subsampling human specific bams
process target_human {
	
	tag "target Hs ${pair_id}"
	label 'big_mem'

	publishDir "${params.outdir}/$pair_id"

	input:
	tuple val(pair_id), path(sorted_dup_bam), path(dup_metrics_txt)
	path refdir

	output:
	tuple val(pair_id), path("${pair_id}.sorted.dup.hs.bam")

	script:
	"""
	samtools view -b -h ${sorted_dup_bam} \
	-T $refdir/genome.fa \
	-L $refdir/human.bed > ${pair_id}.sorted.dup.hs.bam
	"""	
}

// insert size calculation
process insert_sizes {
	
	tag "insert sizes ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"

	input:
	tuple val(pair_id), path(pf_bam)

	output:
	tuple val(pair_id), path("${pair_id}.insert.txt"), path("${pair_id}.insert2.txt"), path("${pair_id}_histo.pdf")
	
	script:
	"""
	gatk CollectInsertSizeMetrics \
	-I ${pf_bam} \
	-O ${pair_id}.insert.txt \
	-H ${pair_id}_histo.pdf \
	-M 0.05

	awk 'FNR>=8 && FNR<=8 {print \$1,\$3,\$4,\$5,\$6,\$7,\$8,\$10,\$11,\$12,\$13,\$14,\$15,\$16,\$17,\$18,\$19,\$20,\$NF="${pair_id}"}' ${pair_id}.insert.txt > ${pair_id}.insert2.txt
	"""
}

// insert summary
process insert_summary {

	tag "insert summary"

	publishDir params.outdir, mode:'copy'

	input:
	path(insert2_collection)

	output:
	path('InsertSize_Final.tsv')

	script:
	"""
	cat $insert2_collection > InsertSize_Final.tsv
	"""
}

// Total bam statistics by sample
process total_bam_stat_per_sample {
	
	tag "Total bam stat ${pair_id}"

	publishDir "${params.outdir}/$pair_id/stat_dir"
	
	input: 
	tuple val(pair_id), path(sorted_dup_bam), path(dup_metrics_txt)

	output:
	tuple val(pair_id), path("${pair_id}_bamstat_total_final.tsv"), path("${pair_id}_bamstat_total.tsv")

	conda 'bioconda::samtools bioconda::datamash'

	script:
	"""
    samtools stats ${pf_bam} | grep ^SN | cut -f 2- | awk -F"\t" '{print \$2}' > ${pair_id}_bamstat_total.tsv
    
    datamash transpose < ${pair_id}_bamstat_total.tsv | awk -F"\t" -v OFS="\t" '{ \$(NF+1) = "${pair_id}"; print }' > ${pair_id}_bamstat_total_final.tsv
    """
}

// Total bam statistic summary
process total_stat_summary {
	
	tag "Total stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	path(bamstat_total) 

	output:
	path('Bam_stats_Total_Final.tsv') 

	script:
	"""
	cat $bamstat_total  | sed '1irow_total_reads_total	filtered_reads_total	sequences_total	is_sorted_total	1st_fragments_total	last_fragments_total	reads_mapped_total	reads_mapped_and_paired_total	reads_unmapped_total	reads_properly_paired_total	reads_paired_total	reads_duplicated_total	reads_MQ0_total	reads_QC_failed_total	non_primary_alignments_total	total_length_total	total_first_fragment_length_total	total_last_fragment_length_total	bases_mapped_total	bases_mapped_(cigar)_total	bases_trimmed_total	bases_duplicated_total	mismatches_total	error_rate_total	average_length_total	average_first_fragment_length_total	average_last_fragment_length_total	maximum_length_total	maximum_first_fragment_length_total	maximum_last_fragment_length_total	average_quality_total	insert_size_average_total	insert_size_standard_deviation_total	inward_oriented_pairs_total	outward_oriented_pairs_total	pairs_with_other_orientation_total	pairs_on_different_chromosomes_total	percentage_of_properly_paired_reads_(%)_total	sample_id' > Bam_stats_total_Final.tsv
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
    """
}

// Pf bam statistic summary
process pf_stat_summary {
	
	tag "Pf stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	path(bamstat_pf) 

	output:
	path('Bam_stats_pf_Final.tsv') 

	script:
	"""
	cat $bamstat_pf  | sed '1irow_total_reads_pf	filtered_reads_pf	sequences_pf	is_sorted_pf	1st_fragments_pf	last_fragments_pf	reads_mapped_pf	reads_mapped_and_paired_pf	reads_unmapped_pf	reads_properly_paired_pf	reads_paired_pf	reads_duplicated_pf	reads_MQ0_pf	reads_QC_failed_pf	non_primary_alignments_pf	total_length_pf	total_first_fragment_length_pf	total_last_fragment_length_pf	bases_mapped_pf	bases_mapped_(cigar)_pf	bases_trimmed_pf	bases_duplicated_pf	mismatches_pf	error_rate_pf	average_length_pf	average_first_fragment_length_pf	average_last_fragment_length_pf	maximum_length_pf	maximum_first_fragment_length_pf	maximum_last_fragment_length_pf	average_quality_pf	insert_size_average_pf	insert_size_standard_deviation_pf	inward_oriented_pairs_pf	outward_oriented_pairs_pf	pairs_with_other_orientation_pf	pairs_on_different_chromosomes_pf	percentage_of_properly_paired_reads_(%)_pf	sample_id' > Bam_stats_pf_Final.tsv
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
    """
}

// Hs bam statistic summary
process hs_stat_summary {
	tag "Hs stat summary"

	publishDir params.outdir, mode:'copy'

	input:
	path(bamstat_hs)

	output:
	path('Bam_stats_hs_Final.tsv')

	script:
	"""
	cat $bamstat_hs | sed '1irow_total_reads_hs	filtered_reads_hs	sequences_hs	is_sorted_hs	1st_fragments_hs	last_fragments_hs	reads_mapped_hs	reads_mapped_and_paired_hs	reads_unmapped_hs	reads_properly_paired_hs	reads_paired_hs	reads_duplicated_hs	reads_MQ0_hs	reads_QC_failed_hs	non_primary_alignments_hs	total_length_hs	total_first_fragment_length_hs	total_last_fragment_length_hs	bases_mapped_hs	bases_mapped_(cigar)_hs	bases_trimmed_hs	bases_duplicated_hs	mismatches_hs	error_rate_hs	average_length_hs	average_first_fragment_length_hs	average_last_fragment_length_hs	maximum_length_hs	maximum_first_fragment_length_hs	maximum_last_fragment_length_hs	average_quality_hs	insert_size_average_hs	insert_size_standard_deviation_hs	inward_oriented_pairs_hs	outward_oriented_pairs_hs	pairs_with_other_orientation_hs	pairs_on_different_chromosomes_hs	percentage_of_properly_paired_reads_(%)_hs	sample_id' > Bam_stats_hs_Final.tsv
	"""
}

// Pf:Hs read ratio calculation
process pf_hs_ratio_calc {
	
	tag "Pf:Hs ratio"

	publishDir params.outdir, mode:'copy'
	
	input: 
	tuple path(bamsum_pf), path(bamsum_hs)

	output:
	file('Ratios_hs_pf_reads.tsv')

	shell:
	'''
	paste !{bamsum_pf} !{bamsum_hs} | awk -v OFS="\t" '!/per/ {print$40,$7,$47}' | awk '{if($2==0) $2=1}1' |  awk '{if($3==0) $3=1}1' | awk -v OFS="\t" '{print$1,$2,$3,($3/$2)}' | sed '1ireads_mapped_pf, reads_mapped_hs, ratio_hs_pf' > Ratios_hs_pf_reads.tsv
	'''

}

// distribution of Pf read depth by chromosome
process pf_read_depth {
	
	tag "read depth Pf chroms ${pair_id}"
	scratch true

	publishDir "${params.outdir}/$pair_id/stat_dir/by_chrom"

	input:
	tuple val(pair_id), path(pf_bam)
	path refdir

	output:
	file("${pair_id}_read_coverage.tsv")

	afterScript "rm -rf TMP"

	script: 
	"""
	mkdir -p TMP

	samtools index -bc ${pf_bam}

	for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14
	    do
	    	gatk --java-options "-Xmx${params.gatk_memory}g -Xms${params.gatk_memory}g" DepthOfCoverage \
			-R "$refdir/Pf3D7.fasta" \
			-O chr"\$i" \
			--output-format TABLE \
			-L Pf3D7_"\$i"_v3 \
			--omit-locus-table true \
			-I ${pf_bam} --tmp-dir TMP
			sed '\${/^Total/d;}' chr"\$i".sample_summary > tmp.chr"\$i".sample_summary
			awk -F"\t" -v OFS="\t" '{ print \$0, \$(NF) = "chr'\$i'" }' tmp.chr"\$i".sample_summary > chr"\$i".sample2_summary
	    done

	cat *.sample2_summary | awk '!/sample_id/ {print \$0}' | sed '1isample_id	total	mean	third_quartile	median	first_quartile	bases_perc_above_15	chromosome' > ${pair_id}_read_coverage.tsv
	"""
}

// concatenate Pf read depth by chromosome summaries
process pf_read_depth_summary {
	
	tag "read depth Pf chroms summary"

	publishDir params.outdir, mode:'copy'

	input:
	path(read_coverage)

	output:
	file("ReadCoverage_final.tsv")

	script: 
	"""
	cat $read_coverage | awk '!/sample_id/ {print \$0}' | sed '1isample_id	total	mean	third_quartile	median	first_quartile	bases_perc_above_15	chromosome' > ReadCoverage_final.tsv
	"""
}

// Rmd run quality report generation
process run_report {
 
	tag "Run quality report"

	publishDir params.outdir, mode:'copy'

	input:
	path(report_files)
	path rscript

	output:
	file('run_quality_report.html')

	afterScript "rm -rf TMP"

	script:
	"""	
	mkdir -p TMP/INSERT_FILES
	cp $report_files TMP
	mv TMP/*.insert.txt TMP/INSERT_FILES/
    Rscript -e 'rmarkdown::render(input = "$rscript", output_dir = getwd(), params = list(directory = "TMP"))'
	"""
}


workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the these reports in your browser --> \nmultiqc = $params.outdir/multiqc_report.html\nrun quality report = $params.outdir/run_quality_report.html\nnextflow summary = $params.outdir/report.html": "Oops .. something went wrong" )
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
    trimmed_reads_ch = trimreads(read_pairs_ch, params.trimadapter)
    
    // fastqc report 
    fastqc_ch = fastqc(trimmed_reads_ch)
    // multiqc report --
	multiqc(fastqc_ch.collect()) 

    // bwa alignment
    sam_ch = bwa_align(trimmed_reads_ch, params.refdir)

	// sam format converter, clean, sort, and mark duplicates
	sam_format_ch = sam_format_converter(sam_ch, params.refdir)
	sam_clean_ch = sam_clean(sam_format_ch, params.refdir)
	sam_sort_ch = sam_sort(sam_clean_ch, params.refdir)
	sam_dup_ch = sam_duplicates(sam_sort_ch, params.refdir)

	// samtools sorting Pf and human reads
	pf_bam_ch = target_pf(sam_dup_ch, params.refdir)
	hs_bam_ch = target_human(sam_dup_ch, params.refdir)

	// distribution of Pf read depth by chromosome -- 
	pf_read_depth_ch = pf_read_depth(pf_bam_ch, params.refdir) 
	pf_read_depth_summary(pf_read_depth_ch.collect())

	// insert size calculation
	inserts_ch = insert_sizes(pf_bam_ch) 
	insert2_ch = inserts_ch.map{T->[T[2]]} // select *.insert2.txt
	// insert summary -- 
	insert_summary(insert2_ch.collect()) 

	// Total bam statistics by sample
	total_bamstat_ch = total_bam_stat_per_sample(sam_dup_ch)
	total_final_bamstat_ch = total_bamstat_ch.map{T->[T[1]]} // select *_bamstat_total_final.tsv
	// Total bam statistic summary
	total_summary_ch = total_stat_summary(total_final_bamstat_ch.collect()) 

	// Pf bam statistics by sample
	pf_bamstat_ch = pf_bam_stat_per_sample(pf_bam_ch)
	pf_final_bamstat_ch = pf_bamstat_ch.map{T->[T[1]]} // select *_bamstat_pf_final.tsv
	// Pf bam statistic summary
	pf_summary_ch = pf_stat_summary(pf_final_bamstat_ch.collect()) 

	// Hs bam statistics by sample
	hs_bamstat_ch = hs_bam_stat_per_sample(hs_bam_ch)
	hs_final_bamstat_ch = hs_bamstat_ch.map{T->[T[1]]} // select *_bamstat_hs_final.tsv
	// Hs bam statistic summary
	hs_summary_ch = hs_stat_summary(hs_final_bamstat_ch.collect())

	// Pf:Hs read ratio calculation -- 
	summary_ch = pf_summary_ch.combine(hs_summary_ch) //combine 
	pf_hs_ratio_calc(summary_ch) 

	// Rmd run quality report generation -- 
	insert1_ch = inserts_ch.map{T->[T[1]]} // select *.insert.txt
	report_files_ch = summary_ch.combine(insert1_ch.collect()) //combine
	run_report(report_files_ch, params.rscript)
}
