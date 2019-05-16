import pysam
import sys
import dnaio

bam_fname = snakemake.input.bam
fastq_fnames = snakemake.input.fastqs

bamfile = pysam.AlignmentFile(bam_fname, "rb")

bam_reads = {read.query_name:read for read in bamfile.fetch(snakemake.wildcards.chrom)}
bam_readnames = set(bam_reads)

bamfile.close()

fastq_reads = {}

for fname in fastq_fnames:
    with dnaio.open(fname) as fh:
        for record in fh:
            name = record.name.split()[0]
            if name in bam_readnames:
                fastq_reads[name] = record

with open(snakemake.output.r1,"w") as ofh1:
    with open(snakemake.output.r2,"w") as ofh2:
        for n in bam_readnames:
            fastq_read = fastq_reads[n]
            bam_read = bam_reads[n]
            ofh1.write(">%s\n%s\n+\n%s\n"  % (fastq_read.name.split()[0], fastq_read.sequence, fastq_read.qualities))
            ofh2.write(">%s\n%s\n+\n%s\n"  % (bam_read.query_name, bam_read.query_sequence, "".join([chr(x + 33) for x in bam_read.query_qualities])))
