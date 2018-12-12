#!/bin/bash
echo "@@@@@ initiating @@@@@"

echo "@@@ Converting SRA to FASTQ @@@"
# SRA to FASTQ
cd ~/path/sratoolkit.2.9.2-centos_linux64/bin
./fasterq-dump -O ~/path/FASTQ -p  SRR5590217
./fasterq-dump -O ~/path/FASTQ -p  SRR5590219
./fasterq-dump -O ~/path/FASTQ -p  SRR5590221
./fasterq-dump -O ~/path/FASTQ -p  SRR5590222

# Quality Control (using fastqc)
cd ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/NMJ_FASTQ
fastqc SRR5590217.fastq SRR5590219.fastq SRR5590221.fastq SRR5590222.fastq

echo "@@@ Trimmomatic @@@"
# Trimmomatic (default parameters)
cd ~/Trimmomatic-0.38
java -jar trimmomatic-0.38.jar SE -phred33 ~/path/FASTQ/SRR5590221.fastq ~/path/FASTQ/221trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 ~/path/FASTQ/SRR5590222.fastq ~/path/FASTQ/222trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 ~/path/FASTQ/SRR5590217.fastq ~/path/FASTQ/217trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 ~/path/FASTQ/SRR5590219.fastq ~/path/FASTQ/219trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Quality Control (again using fastqc)
cd ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/NMJ_FASTQ
fastqc 217trim.fastq 219trim.fastq 221trim.fastq 222trim.fastq

echo "@@@ Creating Hisat2 Index @@@"
# Align using hisat2
cd ~/path
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/GRCm38.p6.genome.fa.gz  
hisat2-build GRCm38.p4.genome.fa.gz NMJidx #build the hisat2 intex from the reference genome

echo "@@@ Aligning 217 @@@"
#gunzip SRR5590217trim.fastq.gz
hisat2 --dta-cufflinks -x  ~/path/NMJidx ~/path/FASTQ/217trim.fastq -S hisat2_SAM_217
#gzip SRR5590217trim.fastq.gz

echo "@@@ Aligning 219 @@@"
#gunzip SRR5590217trim.fastq.gz
hisat2 --dta-cufflinks -x  ~/path/NMJidx ~/path/FASTQ/219trim.fastq -S hisat2_SAM_219
#gzip SRR5590217trim.fastq.gz

echo "@@@ Aligning 221 @@@"
#gunzip SRR5590217trim.fastq.gz
hisat2 --dta-cufflinks -x  ~/path/NMJidx ~/path/FASTQ/221trim.fastq -S hisat2_SAM_221
#gzip SRR5590217trim.fastq.gz

echo "@@@ Aligning 222 @@@"
#gunzip SRR5590219trim.fastq.gz
hisat2 --dta-cufflinks -x ~/path//NMJidx ~/path/FASTQ/222trim.fastq -S hisat2_SAM_222
# hisat2_SAM_alignments: file that alignments will be written to (SAM files)
# --dta-cufflinks formats the alignment for later DE analysis using Cufflinks
#gzip SRR5590219trim.fastq.gz

echo "@@@@@ Checking Alignments @@@@@"
## Check alignments
#cd ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/
samtools quickcheck hisat2_SAM_217.sam
samtools quickcheck hisat2_SAM_219.sam
samtools quickcheck hisat2_SAM_221.sam
samtools quickcheck hisat2_SAM_222.sam

echo "@@@ Converting SAM to BAM files @@@"
# Convert to SAM files to BAM files. This helps reduce amount of storage used for this pipeline if you delete the SAM files
samtools view -bS -o hisat2_alignment_217.bam hisat2_SAM_217
samtools view -bS -o hisat2_alignment_219.bam hisat2_SAM_219
samtools view -bS -o hisat2_alignment_221.bam hisat2_SAM_221
samtools view -bS -o hisat2_alignment_222.bam hisat2_SAM_222

echo "@@@ Sorting BAM Files @@@"
# Sort SAM files on chromosomal alignments
# http://www.htslib.org/doc/samtools.html
samtools sort -l 0 -o hisat2_sorted_217.bam -O bam hisat2_alignment_217.bam
samtools sort -l 0 -o hisat2_sorted_219.bam -O bam hisat2_alignment_219.bam
samtools sort -l 0 -o hisat2_sorted_221.bam -O bam hisat2_alignment_221.bam
samtools sort -l 0 -o hisat2_sorted_222.bam -O bam hisat2_alignment_222.bam

echo "@@@ Retreiving Headers @@@"
# retreive the headers
samtools view -H hisat2_sorted_217.bam

### We did not end up running cufflinks and skipped to HTSeq to generate the count tables for DESeq
echo "@@@ Running Cufflinks @@@"
# run Cufflinks
cufflinks -o 217_clout_2 hisat2_sorted_217.bam
cufflinks -o 219_clout_2 hisat2_sorted_219.bam

echo "@@@ Running Cuffmerge  @@@"
cuffmerge -g gencode.vM19.chr_patch_hapl_scaff.annotation.gtf assemblies.txt

echo "@@@ Running Cuffdiff @@@"
cuffdiff -o diff_out -u merged_asm/merged.gtf -L 217,219 -v  hisat2_sorted_217.bam hisat2_sorted_219.bam

### HTSeq (convert sorted files into count matrix for input into DESeq2
htseq-count -f bam hisat2_sorted_217.bam hisat2_sorted_219.bam gencode.vM19.chr_patch_hapl_scaff.annotation.gtf
