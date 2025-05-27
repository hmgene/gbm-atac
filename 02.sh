bw=(
bigdata/nf/out/bwa/merged_library/bigwig/C16_REP1.mLb.clN.bigWig	
bigdata/nf/out/bwa/merged_library/bigwig/LUT_REP1.mLb.clN.bigWig
bigdata/nf/out/bwa/merged_library/bigwig/DMSO_REP1.mLb.clN.bigWig
)

#mkdir -p data
#hm ucsc-refflat2bed12 <(hm ucsc-refflat hg38 | awk -v OFS="\t" '{$2=$1;print $0;}' ) | cut -f 1-6 | sort -u | tee data/gene.bed |\
#hm bed5p - | sort -u | tee data/tss.bed | hm bedw - 3k > data/tss_3k.bed 
computeMatrix scale-regions -S ${bw[@]} -R data/tss_3k.bed  -o bigdata/tss_3k.matrix.gz 
mkdir -p results
plotHeatmap -m bigdata/tss_3k.matrix.gz -o results/tss_3k.png
