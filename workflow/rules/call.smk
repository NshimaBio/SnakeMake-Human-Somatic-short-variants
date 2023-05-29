rule Mutect:
    input:
        ref=gatk_dict["ref"],
        sam=["results/prepared/{s}.T.cram","results/prepared/{s}.N.cram"],
        sam_idx=["results/prepared/{s}.T.cram.crai","results/prepared/{s}.N.cram.crai"],
        intervals=config["fastq"].get("restrict_regions",""),
        germline=gatk_dict["agnomd"],
        agnomd_idx=gatk_dict["agnomd_idx"],
        pon=gatk_dict["g1k_pon"],
        pon_idx=gatk_dict["g1k_pon_idx"]
    output:
        vcf="results/called/{s}.vcf.gz",
        # sam="results/called/{s}.bam",
        f1r2="results/called/{s}.f1r2.tar.gz",
    threads: 32
    resources:
        mem_mb=2048
    params:
        normal_sample_name="{s}.N"
    log:
        "logs/call/Mutect2/{s}.log"
    # wrapper:
    #     config["warpper_mirror"]+"bio/gatk/mutect"
    conda: "../envs/Mutect.yaml"
    script:
        "../scripts/Mutect.py"

rule GetPileupSummaries:
    input:
        bam="results/prepared/{s}.{g}.cram",
        intervals=config["fastq"].get("restrict_regions",""),
        variants=gatk_dict["small_exac_common_3"],
    output:
        "results/called/{s}.{g}.getpileupsummaries.table"
    threads: 32
    resources:
        mem_mb=1024
    params:
        extra="-R {}".format(gatk_dict["ref"])
    log:
        "logs/call/GetPileupSummaries/{s}.{g}.log"
    wrapper:
        config["warpper_mirror"]+"bio/gatk/getpileupsummaries"

rule CalculateContamination:
    input:
        getpileupsummaries_table="results/called/{s}.T.getpileupsummaries.table",
        matched_normal="results/called/{s}.N.getpileupsummaries.table"
    output:
        contamination="results/called/{s}.contamination.table",
        segmentation="results/called/{s}.segmentation.table",
    conda:
        "../envs/CalculateContamination.yaml"
    params:
        # extra=f"-R {gatk_dict["ref"]}",
    log:
        "logs/call/CalculateContamination/{s}.log",
    script:
        "../scripts/CalculateContamination.py"

rule LearnReadOrientationModel:
    input:
        f1r2="results/called/{s}.f1r2.tar.gz",
    output:
        "results/called/{s}.artifacts_prior.f1r2.tar.gz"
    resources:
        mem_mb=1024
    params:
        extra=""
    log:
        "logs/call/LearnReadOrientationModel/{s}.T.log"
    wrapper:
        config["warpper_mirror"]+"bio/gatk/learnreadorientationmodel"

rule FilterMutectCalls:
    input:
        vcf="results/called/{s}.vcf.gz",
        ref=gatk_dict["ref"],
        intervals=config["fastq"].get("restrict_regions",""),
        contamination="results/called/{s}.contamination.table",
        segmentation="results/called/{s}.segmentation.table",
        f1r2="results/called/{s}.artifacts_prior.f1r2.tar.gz"
    output:
        vcf="results/called/{s}.filtered.vcf.gz"
    log:
        "logs/call/FilterMutectCalls/{s}.log"
    params:
        extra="--max-alt-allele-count 3"  # optional arguments, see GATK docs
    resources:
        mem_mb=1024
    wrapper:
        config["warpper_mirror"]+"bio/gatk/filtermutectcalls"

rule merge_filtered:
    input:
        vcfs=expand("results/called/{s}.filtered.vcf.gz",s=samples.Sample),
    output:
        "results/called/all.vcf.gz",
    log:
        "logs/call/merge_filtered.log",
    resources:
        mem_mb=2048
    wrapper:
        config["warpper_mirror"]+"bio/picard/mergevcfs"
