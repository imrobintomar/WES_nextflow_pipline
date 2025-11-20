process HC {
    tag "${sample_id}:chr${chr}"
    maxForks 24

    input:
    tuple val(sample_id), path(bam)

    // SCATTER CHROMOSOMES — must use "each chr from"
    each chr from ((1..22) + ['X'])

    output:
    tuple val(sample_id), val(chr), path("${sample_id}.${chr}.vcf")

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.ref} \
        -I ${bam} \
        -O ${sample_id}.${chr}.vcf \
        -L chr${chr} \
        --native-pair-hmm-threads ${task.cpus}
    """
}
