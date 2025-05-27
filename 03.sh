

peaks=(
bigdata/nf/out/bwa/merged_library/macs2/broad_peak/C16_REP1.mLb.clN_peaks.broadPeak	
bigdata/nf/out/bwa/merged_library/macs2/broad_peak/LUT_REP1.mLb.clN_peaks.broadPeak
bigdata/nf/out/bwa/merged_library/macs2/broad_peak/DMSO_REP1.mLb.clN_peaks.broadPeak
)

bams=(
bigdata/nf/out/bwa/merged_library/C16_REP1.mLb.clN.sorted.bam	
bigdata/nf/out/bwa/merged_library/LUT_REP1.mLb.clN.sorted.bam
bigdata/nf/out/bwa/merged_library/DMSO_REP1.mLb.clN.sorted.bam
)
fn(){
    cat $1 | awk '$9>-log(0.001)/log(10)' | cut -f1-3 
};export -f fn
#parallel --line-buffer fn {} ::: ${peaks[@]} | sort -k1,1 -k2,3n | bedtools merge -i stdin |\
#bedtools multicov -bams ${bams[@]} -bed stdin > results/merged.counts

edger(){
Rscript -e '
library(data.table)
library(edgeR)
tt=fread("results/merged.counts",col.names=c("chrom","start","end","C16","LUT","DMSO"))

d <- as.matrix(tt[, .(C16, LUT, DMSO)]) 
rownames(d) <- paste0(tt$chrom, ":", tt$start, "-", tt$end)

group <- factor(c("C16", "LUT", "DMSO"))  # Group labels
y <- DGEList(counts=d, group=group)
y <- calcNormFactors(y, method="TMM")  # TMM = Trimmed Mean of M-values

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
dispersion <- 0.1  # or try 0.05 or 0.2 based on prior knowledge

fit <- glmFit(y, design, dispersion=dispersion)
contrast <- makeContrasts(C16vsDMSO = C16 - DMSO, levels=design)
lrt <- glmLRT(fit, contrast=contrast[, "C16vsDMSO"])
y1=lrt$table[lrt$table$PValue < 0.001,]
colnames(y1)=paste0(colnames(y1),"_C16vsDMSO")
y1$id = row.names(y1)

contrast <- makeContrasts(LUTvsDMSO = LUT - DMSO, levels=design)
lrt <- glmLRT(fit, contrast=contrast[, "LUTvsDMSO"])
y2=lrt$table[lrt$table$PValue < 0.001,]
colnames(y2)=paste0(colnames(y2),"_LUTvsDMSO")
y2$id= row.names(y2)

res=merge(y1,y2,all=T)
tmp <- do.call(rbind, strsplit(res$id, "[:-]"))
res$chrom <- tmp[,1]
res$start <- as.integer(tmp[,2])
res$end   <- as.integer(tmp[,3])
res <- res[, c("chrom", "start", "end", setdiff(names(res), c("chrom", "start", "end")))]

res$class <- paste0(
  "C16_", with(res, ifelse(is.na(logFC_C16vsDMSO), "NA", ifelse(logFC_C16vsDMSO > 1, "Up", ifelse(logFC_C16vsDMSO < -1, "Down", "NoChange")))),
  "_LUT_", with(res, ifelse(is.na(logFC_LUTvsDMSO), "NA", ifelse(logFC_LUTvsDMSO > 1, "Up", ifelse(logFC_LUTvsDMSO < -1, "Down", "NoChange"))))
)
fwrite(file="results/merged.p-3.tsv",res,sep="\t")
'
}
#edger
n=results/merged.p-3

#hm cutn $n.tsv chrom,start,end,class > results/tmp
#annotatePeaks.pl results/tmp hg38 -annStats $n.anno.stats  -go results/$n.go  > $n.anno.tsv 

for c in C16_Up_LUT_Up C16_Up_LUT_NA C16_NA_LUT_Up;do
    mkdir -p results/$c
    hm cutn results/merged.p-3.anno.tsv Chr,Start,End,Gene_Name,PeakID | tail -n+2 |\
    grep $c | annotatePeaks.pl - hg38 -go results/$c > results/$c/$c.tsv
done


