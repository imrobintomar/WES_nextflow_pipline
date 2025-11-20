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

    Channel
        .fromFilePairs( params.input, flat: true )
        .map { sample, fq -> 
            def id = sample.tokenize('_')[0]
            tuple(id, fq[0], fq[1])
        }
        .set { READS }

    FASTP( READS )
        .set { CLEAN_READS }

    BWA( CLEAN_READS )
        .set { SAM }

    SORT_MARKDUP( SAM )
        .set { BAM }

    BQSR( BAM )
        .set { RECAL_BAM }

    // scatter 23 chromosomes
    HC( RECAL_BAM )
        .set { CHROM_VCFS }

    // merge chromosome VCFs per sample
    MERGE_VCF( CHROM_VCFS )
        .set { RAW_VCF }

    // SnpSift 1000g annotation (parallel)
    SNP1000G( RAW_VCF )
        .set { ANNOT_1000G }

    // ANNOVAR annotation
    ANNOVAR( ANNOT_1000G )
        .set { ANNOVAR_VCF }

    // final AWK + HEADER filtering
    FINAL_FILTER( ANNOVAR_VCF )
}
