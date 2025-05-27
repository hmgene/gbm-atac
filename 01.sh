
g=( NUP98 RAD51 BRCA1 BRCA2 RAD54L XRCC2 EME1 BLM)
bw=( 
    bigdata/nf/out/bwa/merged_library/bigwig/C16_REP1.mLb.clN.bigWig   
    bigdata/nf/out/bwa/merged_library/bigwig/LUT_REP1.mLb.clN.bigWig
    bigdata/nf/out/bwa/merged_library/bigwig/DMSO_REP1.mLb.clN.bigWig
)

cat `hm home`/data/hg38_gene.bed | grep -wf <( echo ${g[@]} | tr " " "\n" ) |\
tee target.bed  | awk -v OFS="\t" -v w=100000 '{ print $1,$2-w,$3+w,$4,0,$5;}' > target_100k.bed

