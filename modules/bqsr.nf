process BQSR {
    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.recal.bam")

    script:
    """
    gatk BaseRecalibrator \
        -R ${params.ref} \
        -I $bam \
        --known-sites ${params.knownsites} \
        -O ${sample_id}.table

    gatk ApplyBQSR \
        --bqsr-recal-file ${sample_id}.table \
        -I $bam \
        -O ${sample_id}.recal.bam
    """
}
