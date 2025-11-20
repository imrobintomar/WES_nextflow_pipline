process HC {
    tag "$sample_id"
    maxForks 24

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), val(chr), path("${sample_id}.${chr}.vcf")

    each chr in (1..22) + ['X']

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.ref} \
        -I $bam \
        -O ${sample_id}.${chr}.vcf \
        -L chr${chr} \
        --native-pair-hmm-threads ${task.cpus}
    """
}
