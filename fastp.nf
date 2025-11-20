process FASTP {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.fq.gz"), path("${sample_id}_R2.fq.gz")

    script:
    """
    fastp -i $read1 -I $read2 \
        --disable_length_filtering \
        --qualified_quality_phred 30 \
        -o ${sample_id}_R1.fq.gz \
        -O ${sample_id}_R2.fq.gz \
        --html ${sample_id}.html \
        -w ${task.cpus}
    """
}
