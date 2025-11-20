process BWA {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bwa-mem2 mem -t ${task.cpus} -Y -K 100000000 \
        -R '@RG\\tID:0\\tSM:${sample_id}\\tPL:ILLUMINA' \
        ${params.ref} $r1 $r2 > ${sample_id}.sam
    """
}
