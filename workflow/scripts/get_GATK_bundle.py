__author__ = "Victor Wang"
__copyright__ = "Copyright 2023, Victor Wang"
__email__ = "victor@bioquest.cn"
__license__ = "Apache License 2.0"

from snakemake.shell import shell
from pathlib import Path

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
output_dir=Path(snakemake.output[0]).parent
shell(
'''
(gsutil -m cp -r \
    "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz" \
    "gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi" \
    "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz" \
    "gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi" \
    "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz" \
    "gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi" \
    "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict" \
    "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta" \
    "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai" \
    "gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf" \
    {output_dir}) {log}
'''
)

shell(
    """
    (cd {output_dir}
    bgzip -f -@ {snakemake.threads} \
    Homo_sapiens_assembly38.dbsnp138.vcf) {log}
    """
)
# bgzip -f -@ 20 Homo_sapiens_assembly38.dbsnp138.vcf > Homo_sapiens_assembly38.dbsnp138.vcf.gz
# tabix -f -p vcf Homo_sapiens_assembly38.dbsnp138.vcf.gz > Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi
shell(
    """
    (cd {output_dir}
    tabix -f -p vcf Homo_sapiens_assembly38.dbsnp138.vcf.gz) {log}
    """
)

