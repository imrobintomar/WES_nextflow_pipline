process MERGE_VCF {
    input:
    tuple val(sample_id), val(chr), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.merged.vcf")

    script:
    """
    gatk MergeVcfs \
        -I ${sample_id}.chr1.vcf \
        -I ${sample_id}.chr2.vcf \
        -I ${sample_id}.chr3.vcf \
        -I\
        -O ${sample_id}.merged.vcf
    """
}
