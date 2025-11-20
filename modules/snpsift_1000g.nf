process SNP1000G {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.1000g.vcf")

    script:
    """
    java -jar SnpSift.jar annotate \
        /media/drprabudh/m1/annovar/1000/filtered_chrALL.vcf \
        $vcf > ${sample_id}.1000g.vcf
    """
}
