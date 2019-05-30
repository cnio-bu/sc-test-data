import pysam
import sys
import dnaio
import gzip

bam_fname = snakemake.input.bam
fastq_fnames = snakemake.input.fastqs

bamfile = pysam.AlignmentFile(bam_fname, "rb")

bam_reads = {read.query_name:read for read in bamfile.fetch("chr"+snakemake.wildcards.chrom)}
bam_readnames = set(bam_reads)

bamfile.close()

fastq_reads = {}

for fname in fastq_fnames:
    with dnaio.open(fname) as fh:
        for record in fh:
            name = record.name.split()[0]
            if name in bam_readnames:
                fastq_reads[name] = record

with gzip.open(snakemake.output.r1,"w") as ofh1:
    with gzip.open(snakemake.output.r2,"w") as ofh2:
        for n in bam_readnames:
            fastq_read = fastq_reads[n]
            bam_read = bam_reads[n]
            ofh1.write("@{}\n{}\n+\n{}\n".format(fastq_read.name.split()[0], fastq_read.sequence, fastq_read.qualities).encode())
            ofh2.write("@{}\n{}\n+\n{}\n".format(bam_read.query_name, bam_read.query_sequence, "".join([chr(x + 33) for x in bam_read.query_qualities])).encode())
