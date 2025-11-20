/*
 modules/haplotypecaller.nf
 DSL2 module for per-chromosome HaplotypeCaller.
 This module expects input tuples:
   (sample_id, bam_path, chr)
 and emits:
   (sample_id, chr, vcf_path)
*/

process HC {
    tag "${sample_id}:chr${chr}"
    maxForks 24

    // inputs: sample id, bam file, and chromosome identifier (1..22 or 'X')
    input:
    tuple val(sample_id), path(bam), val(chr)

    // output: a tuple containing sample id, chr, and the generated VCF path
    output:
    tuple val(sample_id), val(chr), path("${sample_id}.${chr}.vcf")

    // workflow resources can be set per-process (optional)
    // cpus 4
    // memory '16 GB'

    script:
    """
    set -euo pipefail

    gatk HaplotypeCaller \
        -R ${params.ref} \
        -I ${bam} \
        -O ${sample_id}.${chr}.vcf \
        -L chr${chr} \
        --native-pair-hmm-threads ${task.cpus}
    """
}
