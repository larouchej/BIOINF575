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





















