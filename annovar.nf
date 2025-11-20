process ANNOVAR {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.annovar.vcf")

    script:
    """
    table_annovar.pl $vcf /media/drprabudh/m1/annovar/hg38_humandb \
       --buildver hg38 \
       --out ${sample_id}.annovar \
       --remove \
       --protocol refGeneWithVer,dbnsfp42a,clinvar_20240416,gnomad40_exome,avsnp150,cosmic84_coding,exac03 \
       --operation gx,f,f,f,f,f,f \
       --vcfinput \
       --thread ${task.cpus}
    """
}
