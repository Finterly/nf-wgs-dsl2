/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

params.refdir = "$projectDir/refs/genomes" 
params.rscript = "$projectDir/refs/run_quality_report.Rmd" 
params.reads = "${params.inputdir}/*_R{1,2}*.fastq.gz" // if start from QC 
params.bams = "${params.inputdir}/*.sorted.dup.pf.{bam,bam.csi}" // if start from GVCF

log.info """\
W G S - P I P E L I N E!
================================
qc_only         : $params.qc_only
gvcf_only       : $params.gvcf_only
inputdir        : $params.inputdir
outputdir       : $params.outputdir
trimadapter     : $params.trimadapter
"""

// workflows 
include { QC } from './workflows/qc.nf'
include { GVCF } from './workflows/gvcf.nf'

workflow {
    if(params.qc_only && params.gvcf_only){
        // check parameters
        error "Error: only one of (qc_only, gvcf_only) can be enabled."
    } else if (params.qc_only){
        // qc only
        QC()
    } else if (params.gvcf_only) {
        // gvcf only
        pf_bam_ch = Channel.fromFilePairs(params.bams, checkIfExists: true).map{index, bam_index -> [index, *bam_index.flatten()]}
        GVCF(pf_bam_ch)
    }
    else {
        // qc then gvcf
        QC | GVCF
    }    
}