process FINAL_FILTER {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    path("${sample_id}.filtered.tsv")

    script:
    """
    awk '... (your conditions) ...' $vcf > ${sample_id}.filtered.tsv
    """
}
