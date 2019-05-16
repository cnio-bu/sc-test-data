from snakemake.remote import FTP
from snakemake.remote import HTTP

FTP = FTP.RemoteProvider()
HTTP = HTTP.RemoteProvider()

configfile: "config.yaml"

rule all:
    input:
        expand(["ref/annotation.chr{chrom}.gtf",
                "ref/genome.chr{chrom}.fa"], chrom=config["chrom"]),
        expand("reads/{sample}.chr{chrom}.{group}.fastq", 
               group=[1, 2], sample=config["samples"], chrom=config["chrom"])

rule annotation:
    input:
        FTP.remote(config["ref"]["annotation"], static=True, keep_local=True)
    output:
        "ref/annotation.chr{chrom}.gtf"
    shell:
        "zgrep -P ^{wildcards.chrom} {input} > {output}"

rule genome:
    input:
        FTP.remote(config["ref"]["genome"], static=True, keep_local=True)
    output:
        "ref/genome.chr{chrom}.fa"
    shell:
        "gzip -d -c {input} > {output}"

rule sample_bam:
    output:
        temp("reads/{sample}.chr{chrom}.2.bam")
    params:
        url=lambda wc: config["samples"][wc.sample]["bam"],
        seed=lambda wc: abs(hash(wc.sample)) % 10000
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -b -s{params.seed}{config[sampling_factor]} {params.url} chr{wildcards.chrom} > {output}"

rule gen_fastqs:
    input:
        bam="reads/{sample}.chr{chrom}.2.bam",
        fastqs=lambda wc: HTTP.remote(config["samples"][wc.sample]["fastqs"], static=True, keep_local=True)
    output:
        r1="reads/{sample}.chr{chrom}.1.fastq",
        r2="reads/{sample}.chr{chrom}.2.fastq"
    conda:
        "envs/fastqs.yaml"
    script:
        "scripts/gen_fastqs.py"
