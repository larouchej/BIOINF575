#!/bin/bash

#SBATCH --mem=40g
#SBATCH --cpus-per-task=4
#SBATCH --time=15:00:00

echo "@@@@@ initiating @@@@@"

echo "@@@ Converting SRA to FASTQ @@@"
### SRA to FASTQ
cd ~/path/sratoolkit.2.9.2-centos_linu64/bin
./fasterq-dump -O ~/path/FASTQ -p SRR5590217
./fasterq-dump -O ~/path/FASTQ -p SRR5590219
./fasterq-dump -O ~/path/FASTQ -p SRR5590221
./fasterq-dump -O ~/path/FASTQ -p SRR5590222

### Quality Control (fastqc)
cd <path with FASTQ files>
fastqc SRR5590217.fastq.gz SRR5590219.fastq.gz SRR5590221.fastq.gz SRR5590222.fastq.gz

echo "@@@ Trimmomatic @@@"
### Trimmomatic (default parameters)
cd ~/path/Trimmomatic-0.38
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/217trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/219trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/221trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/222trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#gzip SRR5590217trim.fastq
#gzip SRR5590219trim.fastq
#gzip SRR5590221trim.fastq
#gzip SRR5590222trim.fastq

### Quality Control (fastqc)
#cd ~
#fastqc <path to trimmed fastq files>

### Indexing the reference genome:
STAR --runMode genomeGenerate --genomeDir <path to reference index> GenomeDir --genomeFastaFiles <path of reference sequence FASTA files to be indexed> --sjdbGTFfile <path of GTF files> --sjdbOverhang 100 --runThreadN 4

### Aligning to reference genome:
STAR --genomeDir <path to output files> --readFilesIn <path to fastq files> --outSAMstrandField intronMotif --runThreadN 4

### Convert output SAM to BAM files
#cd <path to SAM files>
samtools view -bS -o <name of output BAM files> <name of input SAM files>

### Sort and indexing aligned reads
#samtools sort 217.bam -o 217.sorted bam
#samtools sort 219.bam -o 219.sorted bam
#samtools sort 221.bam -o 221.sorted bam
#samtools sort 222.bam -o 222.sorted bam

#samtools index 217.sorted.bam
#samtools index 219.sorted.bam
#samtools index 221.sorted.bam
#samtools index 222.sorted.bam

### Counting aligned reads (RStudio, HTSeq)
setwd(<path to sorted.bam files>)
directory = getwd()

htseq-count -f bam <path to sorted.bam files>

### Differential expression analysis (RStudio, DESeq)
BiocManager::install("DESeq", version = "3.8")
BiocManager::install("gplots")

library("DESeq")

### Imports unfiltered counts file and creates count data table
txi.star.counts = read.table(file = '<file names of count tables>', sep ="\t", header=TRUE, row.names=1)
star.countData = sapply(txi.star.counts, as.integer, USE.NAMES=TRUE)
rownames(star.countData) = rownames(txi.star.counts)
head(star.countData)

### Makes samples table
samples = data.frame(row.names = colnames(txi.star.counts),
                     "Condition" = c("Aged", "Adult", "Adult", "Aged"),
                     "Sample" = c("Syn", "Syn", "Ctrl", "Ctrl))
samples$group = factor(paste0(samples$Condition, samples$Sample))
samples

### Create the count dataset
condition = samples$group
cds = newCountDataSet(star.countData, condition)

### Pre-filtering
keep <- rowSums(counts(cds)) >= 1
sum(keep)
lenght(rownames(cds))
fract_kept = sum(keep)/lenght(rownames(cds))
fract_kept
cds <- cds[keep,]

### Perform differential expression
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds, normalized=TRUE))

cds = estimateDispersion(cds, method="blinds", sharingMode="fit-only")
str(fitINfo(cds))
plotDispEsts(cds)

res = nbinomTest(cds, "AgedSyn", "AdultSyn")

head(res)
head(res[order(res$padj),])
addmargins(table(res_sig=res$padj<.05))

geneid = (res$id)
geneid = unlist(lapply(geneid, function(x) {
    x <- strsplit(x,split=".",fixed=TRUE)[[1]][1]
    return(x)
}))
length(unique(geneid))==length(rownames(cds))
resTable = (res)
resTable$truncID = geneid

library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)
resTable$symbol = mapIds(org.Mm.eg.db,
                         keys = geneid,
                         column="SYMBOL",
                         keytype="ENSEMBL",
                         multiVals="first")

plotMA(res, ylim=c(-3,3))

ache = match('Ache', resTable$symbol)
musk = match('Musk', resTable$symbol)
hdac4 = match('Hdac4', resTable$symbol)
lrp4 = match('Lrp4', resTable$symbol)
chrna1 = match('Chrna1', resTable$symbol)
myog = match('Myog', resTable$symbol)
myod = match('Myod', resTable$symbol)
ppargc1a = match('Ppargc1a', resTable$symbol)
esrra = match('Esrra', resTable$symbol)
mfn2 = match('Mfn2', resTable$symbol)
opa1 = match('Opa1', resTable$symbol)
pgp = match('Pgp', resTable$symbol)
itgb1 = match('Itgb1', resTable$symbol)
gadd45a = match('Gadd45a', resTable$symbol)
mrpl15 = match('Mrpl15', resTable$symbol)
myl4 = match('Myl4', resTable$symbol)
fst = match('Fst', resTable$symbol)

keep20 <- (ache,musk,hdac4,lrp4,chrna1,myog,myod,ppargc1a,esrra,mfn2,opa1,pgp,itgb1,gadd45a,mrp15,myl4,fst)
keep20
new_cds <- cds[keep20,]

vsd = varianceStabilizingTransformation(new_cds)
library("RColorBrewer")
librar("gplots")
select = order(rowMeans(counts(new_cds)),decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd), col=hmcol, trace="none", margins=c(13,13))





















