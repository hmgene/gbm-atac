wd=bigdata/nf
mkdir -p $wd
cd $wd

echo \
'sample,fastq_1,fastq_2,replicate
C16,../C16_1.fq.gz,../C16_2.fq.gz,1
DMSO,../DMSO_1.fq.gz,../DMSO_2.fq.gz,1
LUT,../LUT_1.fq.gz ,../LUT_2.fq.gz,1' > input.csv

echo "#!/bin/bash
module load singularity
nextflow run nf-core/atacseq -r 2.1.2 --input input.csv -profile singularity --outdir out \
--genome hg38 --read_length 150 --macs_gsize 2700000000
" | sbatch --mem=128g -c 24 --time=100:00:00

cd -
