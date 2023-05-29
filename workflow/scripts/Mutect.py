"""Snakemake wrapper for GATK4 Mutect2"""

__author__ = "Victor Wang"
__copyright__ = "Copyright 2023, Victor Wang"
__email__ = "victor@bioquest.cn"
__license__ = "Apache License 2.0"

import os
import tempfile
from snakemake.shell import shell
from snakemake.utils import makedirs
from snakemake_wrapper_utils.java import get_java_opts

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# On non-omp systems, and in case OMP_NUM_THREADS
# was not set, define OMP_NUM_THREADS through python
if "OMP_NUM_THREADS" not in os.environ.keys() and snakemake.params.get(
    "use_omp", False
):
    os.environ["OMP_NUM_THREADS"] = snakemake.threads


sam_output = snakemake.output.get("sam", "")
if sam_output:
    sam_output = f"--bam-output {sam_output }"


germline_resource = snakemake.input.get("germline", "")
if germline_resource:
    germline_resource = f"--germline-resource {germline_resource}"

intervals = snakemake.input.get("intervals", "")
if intervals:
    intervals = f"--intervals {intervals}"

f1r2 = snakemake.output.get("f1r2", "")
if f1r2:
    f1r2 = f"--f1r2-tar-gz {f1r2}"

pon = snakemake.input.get("pon", "")
if pon:
    pon = f"--panel-of-normals {pon}"

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

# In case Java execution environment suits GC parallel
# calls, these must be given as optional java parameters
if snakemake.params.get("use_parallelgc", False):
    if "UseParallelGC" not in java_opts:
        java_opts += " -XX:+UseParallelGC "
    if "ParallelGCThreads" not in java_opts:
        java_opts += f" -XX:ParallelGCThreads={snakemake.threads}"

sam = snakemake.input.maps
if isinstance(sam,str):
    sam=[sam]
sam = " ".join("--input {}".format(x) for x in sam)

normal = snakemake.params.get("normal_sample_name", "")
if isinstance(normal,str):
    normal=[normal]
normal = " ".join("--normal-sample {}".format(x) for x in normal)

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk --java-options '{java_opts}' Mutect2"  # Tool and its subprocess
        " --native-pair-hmm-threads {snakemake.threads}"
        " {sam}"  # Path to input mapping file
        " --reference {snakemake.input.ref}"  # Path to reference fasta file
        " {normal}"
        " {f1r2}"  # Optional path to output f1r2 count file
        " {germline_resource}"  # Optional path to optional germline resource VCF
        " {intervals}"  # Optional path to optional bed intervals
        " {pon} "  # Optional path to panel of normals
        " {extra}"  # Extra parameters
        " --tmp-dir {tmpdir}"
        " --output {snakemake.output.vcf}"  # Path to output vcf file
        " {bam_output}"  # Path to output bam file, optional
        " {log}"  # Logging behaviour
    )
