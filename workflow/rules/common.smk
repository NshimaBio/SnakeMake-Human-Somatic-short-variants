import numpy as np
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("7.25.0")

report: "../report/workflow.rst"

# container: "mambaorg/micromamba:1.4.2"

#=====================================================
# validate config.yaml file and samples.csv file
#=====================================================

configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")

# import yaml
# with open('config/config.yaml', 'r') as file:
#     config = yaml.safe_load(file)
samples = pd.read_csv(config["samples"], dtype=str,sep='\t',header=0).fillna(value="")
if not "Unit" in samples.columns:
    samples.loc[:,"Unit"]=""
else:
    samples.loc[:,"Unit"] = [f".{x}" if x else "" for x in samples.Unit]

samples.set_index(keys=["Sample", "Unit","Group"], drop=False,inplace=True)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])

validate(samples, schema="../schemas/samples.schema.yaml")

dir = config["fastq"].get("dir")
if config["fastq"].get("pe"):
    fq1=[f"{dir}/{x}{y}.{z}.1.fastq.gz" for x,y,z in zip(samples.Sample,samples.Unit,samples.Group)]
    fq2=[f"{dir}/{x}{y}.2.fastq.gz" for x,y,z in zip(samples.Sample,samples.Unit,samples.Group)]
    samples.insert(loc=0,column="fq2",value=fq2)
    samples.insert(loc=0,column="fq1",value=fq1)
else:
    fq1=[f"{dir}/{x}{y}.{z}.fastq.gz" for x,y,z in zip(samples.Sample,samples.Unit,samples.Group)]
    samples.insert(loc=0,column="fq1",value=fq1)

#=====================================================
# Helper functions
#=====================================================
def get_vep_prefix(config=config):
    g=config["GATK"]
    p="{0}{1}_{2}_{3}_".format(g['dir'],g['species'].capitalize(),g['build'],g['release'])
    return p
def get_dedup_extra(config=config):
    if config["fastq"]["pe"]:
        return "-M"
    else:
        return "-M --ignoreUnmated"

# def get_intervals(config=config,defalut=""):
#     rr=config["fastq"]["restrict_regions"]
#     if rr:
#         rr=pd.read_csv(rr,sep='\t',header=None,dtype=str)
#     else:
#         return defalut

def get_fastq(wildcards):
	"""Get fastq files of given sample and unit."""
	fastqs = samples.loc[(wildcards.s, wildcards.u, wildcards.g), ]
	if config["fastq"].get("pe"):
		return [fastqs.fq1, fastqs.fq2]
	return [fastqs.fq1]

def get_trimmed_fastq(wildcards):
    """Get trimmed reads of given sample and unit."""
    if config["fastq"].get("pe"):
        # paired-end sample
        return expand("results/trimmed/{{s}}{{u}}.{{g}}.{_}.fastq.gz",_=["1","2"])
    # single end sample
    return "results/trimmed/{s}{u}.{g}.fastq.gz"

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{s}{u}.{g}\tSM:{s}.{g}\tLB:{s}.{g}\tPL:{pl}'".format(
        s=wildcards.s,
        u=wildcards.u,
        g=wildcards.g,
        pl=samples.loc[(wildcards.s,wildcards.u,wildcards.g),"Platform"]
    )

def get_ApplyBQSR_cram(wildcards):
    """Get all aligned reads of given sample."""
    u=samples.loc[wildcards.s,"Unit"]
    if u.all(axis=None):
        return expand(
            ["results/prepared/{{s}}{u}.{{g}}.ApplyBQSR.cram"],
            # s=wildcards.s,
            u=u
        )
    else:
        return "results/prepared/{s}.{g}.ApplyBQSR.cram"

# samples.swaplevel(axis=0).loc["P1","Unit"]
def get_cram(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/prepared/{s}.{g}.cram",
        s=wildcards.s,
        g=samples.loc[wildcards.s,"Group"]
    )

def get_cram_idx(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/prepared/{s}.{g}.cram.crai",
        s=wildcards.s,
        g=samples.loc[wildcards.s,"Group"]
    )

# def get_normal_sample_name(wildcards,default=""):
#     """Get the matched normal sample name of given sample."""
#     if wildcards.g:
#         return expand("--normal-sample {s}.N",
#             s=wildcards.s,
#             u=samples.loc[wildcards.s,"Unit"]
#         )
#     else:
        # return default

# def get_interval_padding(config=config,default=""):
#     padding = config["HaplotypeCaller"].get("interval_padding","")
#     if padding:
#         return "--interval-padding {}".format(padding)
#     return default

gatk_dict=dict(ref=config["GATK"]["dir"]+"Homo_sapiens_assembly38.fasta",
        fai=config["GATK"]["dir"]+"Homo_sapiens_assembly38.fasta.fai",
        dict=config["GATK"]["dir"]+"Homo_sapiens_assembly38.dict",
        g1k_pon=config["GATK"]["dir"]+"1000g_pon.hg38.vcf.gz",
        g1k_pon_idx=config["GATK"]["dir"]+"1000g_pon.hg38.vcf.gz.tbi",
        agnomd=config["GATK"]["dir"]+"af-only-gnomad.hg38.vcf.gz",
        agnomd_idx=config["GATK"]["dir"]+"af-only-gnomad.hg38.vcf.gz.tbi",
        small_exac_common_3=config["GATK"]["dir"]+"small_exac_common_3.hg38.vcf.gz",
        small_exac_common_3_idx=config["GATK"]["dir"]+"small_exac_common_3.hg38.vcf.gz.tbi",
        dbsnp=config["GATK"]["dir"]+"Homo_sapiens_assembly38.dbsnp138.vcf.gz",
        dbsnp_idx=config["GATK"]["dir"]+"Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
)

def get_GTAK_bundle_output():
    return multiext(config["GATK"]["dir"],
            "Homo_sapiens_assembly38.fasta",
            "Homo_sapiens_assembly38.fasta.fai",
            "Homo_sapiens_assembly38.dict",
            "1000g_pon.hg38.vcf.gz",
            "1000g_pon.hg38.vcf.gz.tbi",
            "af-only-gnomad.hg38.vcf.gz",
            "af-only-gnomad.hg38.vcf.gz.tbi",
            "small_exac_common_3.hg38.vcf.gz",
            "small_exac_common_3.hg38.vcf.gz.tbi",
            "Homo_sapiens_assembly38.dbsnp138.vcf.gz",
            "Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
            )
#=====================================================
# Wildcard constraints
#=====================================================

wildcard_constraints:
    v="SNP|INDEL",
    s="|".join(samples.index.get_level_values(0)),
    u="|".join(samples.index.get_level_values(1)),
    g="|".join(samples.index.get_level_values(2))