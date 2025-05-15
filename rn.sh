echo \
'sample,fastq_1,fastq_2,replicate
C16,bigdata/C16_1.fq.gz,bigdata/C16_2.fq.gz,1
DMSO,bigdata/DMSO_1.fq.gz,bigdata/DMSO_2.fq.gz,1
LUT, bigdata/LUT_1.fq.gz ,bigdata/LUT_2.fq.gz,1' > input.csv



nextflow run nf-core/atacseq -r 2.1.2 --input input.csv -profile docker --outdir out \
--genome hg38 --read_length 150 --macs_gsize 2700000000
