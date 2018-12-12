#!/bin/bash
echo "@@@@@ initiating @@@@@"


#echo "@@@ Converting SRA to FASTQ @@@"
# Downloading SRR files
# prefetch <SRR#>
# SRA to FASTQ
#cd ~/path/sratoolkit.2.9.2-centos_linux64/bin
#./fasterq-dump -O ~/path/FASTQ -p  SRR5590217
#./fasterq-dump -O ~/path/FASTQ -p  SRR5590219
#./fasterq-dump -O ~/path/FASTQ -p  SRR5590221
#./fasterq-dump -O ~/path/FASTQ -p  SRR5590222

# Quality Control (using fastqc)
#cd <path to fastq files>
#fastqc SRR5590217.fastq.gz SRR5590219.fastq.gz SRR5590221.fastq.gz SRR5590222.fastq.gz

echo "@@@ Trimmomatic @@@"
# Trimmomatic (default parameters)
cd ~/Trimmomatic-0.38
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/217trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/219trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/221trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 <path to fastq files> ~/path/FASTQ/222trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#gzip SRR5590217trim.fastq
#gzip SRR5590219trim.fastq
#gzip SRR5590221trim.fastq
#gzip SRR5590222trim.fastq

# Quality Control (again using fastqc)
#cd ~
#fastqc <path to trimmed fastq files> *fastq.gz

# Indexing the refrence genome:
# cd <path to where Reference Genome fasta file is located and indexed one will be saved>
bowtie2-build -f GRCm38.p4.genome.fa GRCm38_M6


# aligning to refernce genome:
bowtie2 -q --phred33 -p 8 --no-unal -x < Path to indexed refrence genome GRCm38_M6 > -U <Path to trimmed fastq files > | samtools view -bhS -o out.bam
# run the same code for all other samples

# Sort and indexing aligned reads
#samtools sort 217.bam -o 217.sorted.bam
#samtools sort 219.bam -o 219.sorted.bam
#samtools sort 221.bam -o 221.sorted.bam
#samtools sort 222.bam -o 222.sorted.bam

#samtools index 217.sorted.bam
#samtools index 219.sorted.bam
#samtools index 221.sorted.bam
#samtools index 222.sorted.bam

# Counting aligned reads and making a count matrix using Rsubread in R
R

library("Rsubread")
library("DESeq2")
setwd("path to sorted bam files")
filenames = dir(path='path to sorted bam files',
                pattern = "*\\.sorted.bam$",
                full.names=TRUE)

tmp <- Rsubread::featureCounts(filenames, annot.ext="<path to annotation GTF files>",
                     isGTFAnnotationFile = T,
                     ignoreDup=F)

save(tmp, file="./bowtie2.countMatrix.RData")


# Running DEseq on the count matrix and making table for samples
library(DEseq)

load("./bowtie2.countMatrix.RData")

# making the sample table and defining condition for our design

coldata = data.frame("Sample" = c("Syn","Syn","Ctrl","Ctrl"),
                     "Condition" = c("Aged", "Adult", "Adult", "Aged"))
rownames(coldata) = c("AgedSyn", "AdultSyn", "AdultCtrl", "AgedCtrl")
coldata$SRR = c("SRR5590217", "SRR5590219", "SRR5590221", "SRR5590222")
coldata$group = factor(paste0(coldata$Condition, coldata$Sample))
condition = coldata$group
coldata$group

cts <- as.matrix(tmp$counts)
colnames(cts) = c("AgedSyn", "AdultSyn", "AdultCtrl", "AgedCtrl")

cds = newCountDataSet(cts, condition )

keep <- rowSums(counts(cds)) >= 1
sum(keep)
length(rownames(cds))
fract_kept = sum(keep)/length(rownames(cds))
fract_kept
cds <- cds[keep,]
# estimating the size factor and dispersions and plotting the dispersion plot
cds = estimateSizeFactors( cds )
sizeFactors(cds)
head( counts( cds, normalized=TRUE ) )

cds = estimateDispersions(cds, method="blind", sharingMode="fit-only" )
str( fitInfo(cds) )

plotDispEsts( cds )

# running negetive binomial distribution test on the sample, label padj value < .05
res = nbinomTest( cds, "AgedSyn", "AdultSyn" )

head(res)
head(res[order(res$padj),])
addmargins( table( res_sig = res$padj < .05) )
write.csv(res,file="res_correct.csv")

# Gene enseble ids are actually isoform ids and have decimals which make it hard to find their gene counterpart using annotation db, so first I removed the decimals

geneid = res$id
geneid = unlist(lapply(geneid, function(x) {
  x <- strsplit(x,split='.',fixed=TRUE)[[1]][1]
  return(x)
}))
length(unique(geneid))==length(rownames(res)) # making sure that we're not losing any genes
res$id = geneid
head(rownames(res))

# adding gene symbol and enterz id to our result table
ibrary("AnnotationDbi")
library("org.Mm.eg.db")
columns(org.Mm.eg.db)
res$symbol = mapIds(org.Mm.eg.db,
                        keys=res$id,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                         keys=res$id,
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")


# making the heatmap of 20 genes that was in the paper
First_20 <- res[order(res$padj),]
#____ simlar plot as the paper
ache = match('Ache', res$symbol)
musk = match('Musk', res$symbol)
hdac4 = match('Hdac4', res$symbol)
lrp4 = match('Lrp4', res$symbol)
chrna1 = match('Chrna1', res$symbol)
myog = match('Myog', res$symbol)
myod = match('Myod1', res$symbol)
ppargc1a = match('Ppargc1a', res$symbol)
esrra = match('Esrra', res$symbol)
mfn2 = match('Mfn2', res$symbol)
opa1 = match('Opa1', res$symbol)
pgp = match('Pgp', res$symbol)
itgb1 = match('Itgb1', res$symbol)
gadd45a = match('Gadd45a', res$symbol)
mrpl15 = match('Mrpl15', res$symbol)
myl4 = match('Myl4', res$symbol)
fst = match('Fst', res$symbol)

keep20 <- c(ache,musk,hdac4,lrp4,chrna1,myog,myod,ppargc1a,esrra,mfn2,opa1,pgp,itgb1,gadd45a,mrpl15,myl4,fst)
keep20
cds_selected <- cds[keep20,]

# log transforming the data so the mean is not as dependent on variance, also this is necessary to have heatmap
vsd_20 = varianceStabilizingTransformation( cds_selected )
library("RColorBrewer")
library("gplots")
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd_20), col = hmcol, trace="none", margins = c(8,13), Colv = NA, Rowv = NA)

# MA plot, PCA plot and heatmap plot for all genes:
plotMA(res)

vsd = varianceStabilizingTransformation( cds )
library("RColorBrewer")
library("gplots")
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd)[select,], col = hmcol, trace="none", margin=c(13, 13), Colv = FALSE)

print(plotPCA(vsd, intgroup=c("condition")))


# exit R
quit()