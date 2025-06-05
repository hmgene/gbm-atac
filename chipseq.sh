input=(
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM318_DMSO_NKRF_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM318_Luteolin_NKRF_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM318_PKRi_NKRF_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_DMSO_H3K27ac_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_DMSO_NKRF_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_DMSO_p65_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_DMSO_RPB1_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_Luteolin_NKRF_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_Luteolin_p65_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_PKRi_NKRF_R1.fastq.gz
/mnt/vstor/SOM_GENE_BEG33/data/GBM_Maria/GBM358_PKRi_p65_R1.fastq.gz
)

run-nf(){
	odir=bigdata/chipseq; mkdir -p $odir
	cd $odir
	echo sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate  > input.csv
	for i in ${input[@]};do
		s=${i##*/};s=${s%_R1.fastq.gz};
		d=`echo $s| cut -d "_" -f 2`; 
		a=`echo $s| cut -d "_" -f 3`;
		c="$s"; ## itself
		if [ "$d" != "DMSO" ];then 
			c="${s/$d/DMSO}";
		fi
		echo "$s,$i,${i/_R1/_R2},1,$a,$c,1" >> input.csv
	done

	echo "#!/bin/bash
	module load singularity
	nextflow run nf-core/chipseq --input input.csv --outdir ./ -profile singularity \
	--genome hg38 --read_length 150 --macs_gsize 2700000000

	" | sbatch --mem=128g -c 24 --time=100:00:00
	cd -
}


peaks=(
bigdata/chipseq/bwa/merged_library/macs3/GBM318_DMSO_NKRF_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM318_Luteolin_NKRF_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM318_PKRi_NKRF_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
#bigdata/chipseq/bwa/merged_library/macs3/GBM358_DMSO_H3K27ac_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_DMSO_NKRF_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_DMSO_p65_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_DMSO_RPB1_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_Luteolin_NKRF_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_Luteolin_p65_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_PKRi_NKRF_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
bigdata/chipseq/bwa/merged_library/macs3/GBM358_PKRi_p65_REP1.mLb.clN.sorted.bam_peaks.narrowPeak
)

bams=(
bigdata/chipseq/bwa/merged_library/GBM318_DMSO_NKRF_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM318_Luteolin_NKRF_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM318_PKRi_NKRF_REP1.mLb.clN.sorted.bam
#bigdata/chipseq/bwa/merged_library/GBM358_DMSO_H3K27ac_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_DMSO_NKRF_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_DMSO_p65_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_DMSO_RPB1_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_Luteolin_NKRF_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_Luteolin_p65_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_PKRi_NKRF_REP1.mLb.clN.sorted.bam
bigdata/chipseq/bwa/merged_library/GBM358_PKRi_p65_REP1.mLb.clN.sorted.bam
)
o=results/merged_chipseq.counts
printf "chrom\tstart\tend" > $o
echo ${bams[@]} | perl -ne 'while($_=~/(GBM\w+)_REP1/g){ print "\t$1";}'  >> $o
fn(){ cat $1 | awk '$9>=-log(0.001)/log(10)' | cut -f1-3; };export -f fn;
parallel --line-buffer fn {} ::: ${peaks[@]} | sort -k1,1 -k2,3n | bedtools merge -i stdin |\
bedtools multicov -bams ${bams[@]} -bed stdin >> $o


exit


#WT_BCATENIN_IP,BLA203A1_S27_L006_R1_001.fastq.gz,,1,BCATENIN,WT_INPUT,1
#WT_BCATENIN_IP,BLA203A25_S16_L001_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
#WT_BCATENIN_IP,BLA203A25_S16_L002_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
#WT_BCATENIN_IP,BLA203A25_S16_L003_R1_001.fastq.gz,,2,BCATENIN,WT_INPUT,2
#WT_BCATENIN_IP,BLA203A49_S40_L001_R1_001.fastq.gz,,3,BCATENIN,WT_INPUT,3
#WT_INPUT,BLA203A6_S32_L006_R1_001.fastq.gz,,1,,,
#WT_INPUT,BLA203A30_S21_L001_R1_001.fastq.gz,,2,,,
#WT_INPUT,BLA203A30_S21_L002_R1_001.fastq.gz,,2,,,
#WT_INPUT,BLA203A31_S21_L003_R1_001.fastq.gz,,3,,,
