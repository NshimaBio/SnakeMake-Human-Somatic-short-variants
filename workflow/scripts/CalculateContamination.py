__author__ = "Victor Wang"
__copyright__ = "Copyright 2023, Victor Wang"
__email__ = "victor@bioquest.cn"
__license__ = "Apache License 2.0"

import tempfile
from snakemake.shell import shell
from snakemake_wrapper_utils.java import get_java_opts

extra = snakemake.params.get("extra", "")
java_opts = get_java_opts(snakemake)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

matched_normal = snakemake.input.get("matched_normal", "")
if matched_normal:
    matched_normal = f"--matched-normal {matched_normal}"

segmentation = snakemake.output.get("segmentation", "")
if segmentation:
    segmentation = f"--tumor-segmentation {segmentation}"

with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk --java-options '{java_opts}' CalculateContamination"
        " --input {snakemake.input.getpileupsummaries_table}"
        " {matched_normal}"
        " {extra}"
        " --tmp-dir {tmpdir}"
        " --output {snakemake.output.contamination}"
        " {segmentation}"
        " {log}"
    )