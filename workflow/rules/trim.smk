if config["fastq"].get("pe"):
    rule fastp_pe:
        input:
            sample=get_fastq
        output:
            trimmed=[temp("results/trimmed/{s}{u}.{g}.1.fastq.gz"), temp("results/trimmed/{s}{u}.{g}.2.fastq.gz")],
            html=temp("report/{s}{u}.{g}.fastp.html"),
            json=temp("report/{s}{u}.{g}.fastp.json"),
        log:
            "logs/trim/{s}{u}.{g}.log"
        threads: 32
        wrapper:
            config["warpper_mirror"]+"bio/fastp"
else:
    rule fastp_se:
        input:
            sample=get_fastq
        output:
            trimmed="results/trimmed/{s}{u}.{g}.fastq.gz",
            html="report/{s}{u}.{g}.fastp.html",
            json="report/{s}{u}.{g}.fastp.json",
        log:
            "logs/trim/{s}{u}.{g}.log"
        threads: 32
        wrapper:
            config["warpper_mirror"]+"bio/fastp"
