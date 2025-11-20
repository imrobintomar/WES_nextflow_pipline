nextflow.enable.dsl=2

include { FASTP } from './modules/fastp.nf'
include { BWA } from './modules/bwa.nf'
include { SORT_MARKDUP } from './modules/sort_markdup.nf'
include { BQSR } from './modules/bqsr.nf'
include { HC } from './modules/haplotypecaller.nf'
include { MERGE_VCF } from './modules/merge_vcf.nf'
include { SNP1000G } from './modules/snpsift_1000g.nf'
include { ANNOVAR } from './modules/annovar.nf'
include { FINAL_FILTER } from './modules/final_filter.nf'

workflow {

    /*
     * Input FASTQ -> channel of tuples:
     * (sample_id, read1_path, read2_path)
     */
    Channel
        .fromFilePairs( params.input, flat: true )
        .map { sample, fq ->
            def id = sample.tokenize('_')[0]
            tuple(id, fq[0], fq[1])
        }
        .set { READS }

    /*
     * Core preprocessing / alignment chain:
     * FASTP -> BWA -> SORT+MARKDUP -> BQSR
     *
     * Each process is a DSL2 module that accepts and emits tuples:
     * FASTP: (sample_id, r1, r2) -> (sample_id, r1_trim, r2_trim)
     * BWA:   (sample_id, r1_trim, r2_trim) -> (sample_id, sam)
     * SORT_MARKDUP: (sample_id, sam) -> (sample_id, markdup.bam)
     * BQSR:  (sample_id, markdup.bam) -> (sample_id, recal.bam)
     */

    FASTP( READS )
        .set { CLEAN_READS }

    BWA( CLEAN_READS )
        .set { SAM }

    SORT_MARKDUP( SAM )
        .set { BAM }

    BQSR( BAM )
        .set { RECAL_BAM }   // RECAL_BAM is a channel of tuples: (sample_id, recal.bam)

    /*
     * Chromosome channel: produce values 1..22 and X
     * This channel drives the scatter (per-chromosome HaplotypeCaller)
     */
    ch_chrs = Channel.from( (1..22) + ['X'] )

    /*
     * Combine every sample BAM with each chromosome.
     * This produces for each (sample, bam) and each chr -> (sample_id, bam, chr)
     */
    RECAL_BAM
        .combine(ch_chrs)
        .map { sample_tuple, chr ->
            // sample_tuple is (sample_id, bam_path)
            def (sid, bam) = sample_tuple
            tuple(sid, bam, chr)
        }
        .set { BAM_WITH_CHR }

    /*
     * Run HaplotypeCaller per (sample,chr)
     * HC accepts tuples (sample_id, bam, chr) and emits (sample_id, chr, vcf)
     */
    HC( BAM_WITH_CHR )
        .set { CHROM_VCFS }

    /*
     * MERGE the per-chromosome VCFs back to per-sample VCFs
     */
    MERGE_VCF( CHROM_VCFS )
        .set { RAW_VCF }

    /*
     * Annotation and final filtering
     */
    SNP1000G( RAW_VCF )
        .set { ANNOT_1000G }

    ANNOVAR( ANNOT_1000G )
        .set { ANNOVAR_VCF }

    FINAL_FILTER( ANNOVAR_VCF )

} // workflow
