process SORT_MARKDUP {
    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam")

    script:
    """
    gatk SortSam -I $sam -O ${sample_id}.sorted.bam -SO coordinate

    gatk MarkDuplicates \
        -I ${sample_id}.sorted.bam \
        -O ${sample_id}.markdup.bam \
        -M ${sample_id}.dup.txt
    """
}
