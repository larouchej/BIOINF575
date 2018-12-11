#!/bin/bash
echo "@@@@@ initiating @@@@@"

#echo "@@@ Converting SRA to FASTQ @@@"
# SRA to FASTQ
#cd ~/path/sratoolkit.2.9.2-centos_linux64/bin
#./fasterq-dump -O ~/path/FASTQ -p  SRR5590221
#./fasterq-dump -O ~/path/FASTQ -p  SRR5590222

# Quality Control (using fastqc)
#cd ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/NMJ_FASTQ
#fastqc SRR5590217.fastq.gz SRR5590219.fastq.gz

echo "@@@ Trimmomatic @@@"
# Trimmomatic (default parameters)
cd ~/Trimmomatic-0.38
java -jar trimmomatic-0.38.jar SE -phred33 ~/path/FASTQ/SRR5590221.fastq ~/path/FASTQ/221trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar trimmomatic-0.38.jar SE -phred33 ~/path/FASTQ/SRR5590222.fastq ~/path/FASTQ/222trim.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#gzip SRR5590217trim.fastq
#gzip SRR5590219trim.fastq

# Quality Control (again using fastqc)
#cd ~
#fastqc ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/SRR5590217trim.fastq.gz 
#fastqc ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/SRR5590219trim.fastq.gz

#echo "Creating Hisat2 Index @@@"
# Align using hisat2
#cd ~/path
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M19/GRCm38.p6.genome.fa.gz  
#hisat2-build GRCm38.p4.genome.fa.gz NMJidx #build the hisat2 intex from the reference genome

echo "@@@ Aligning 221 @@@"
#gunzip SRR5590217trim.fastq.gz
hisat2 --dta-cufflinks -x  ~/path/NMJidx ~/path/FASTQ/221trim.fastq -S hisat2_SAM_221
#gzip SRR5590217trim.fastq.gz

#129618715 reads; of these:
#  129618715 (100.00%) were unpaired; of these:
#    3441871 (2.66%) aligned 0 times
#    110218282 (85.03%) aligned exactly 1 time
#    15958562 (12.31%) aligned >1 times
#97.34% overall alignment rate

echo "@@@ Aligning 222 @@@"
#gunzip SRR5590219trim.fastq.gz
hisat2 --dta-cufflinks -x ~/path//NMJidx ~/path/FASTQ/222trim.fastq -S hisat2_SAM_222
# hisat2_SAM_alignments: file that alignments will be written to (SAM files)
# --dta-cufflinks formats the alignment for later DE analysis using Cufflinks
#gzip SRR5590219trim.fastq.gz

#94458039 reads; of these:
#  94458039 (100.00%) were unpaired; of these:
#    557396 (0.59%) aligned 0 times
#    85082244 (90.07%) aligned exactly 1 time
#    8818399 (9.34%) aligned >1 times
#99.41% overall alignment rate

echo "@@@@@ Checking Alignments @@@@@"
## Check alignments
#cd ~/Documents/UMichigan/Fall2018/BIOINF575/NMJProject/
samtools quickcheck hisat2_SAM_221.sam
samtools quickcheck hisat2_SAM_222.sam

echo "@@@ Converting SAM to BAM files @@@"
# Convert to SAM files to BAM files
samtools view -bS -o hisat2_alignment_221.bam hisat2_SAM_221
samtools view -bS -o hisat2_alignment_222.bam hisat2_SAM_222

echo "@@@ Sorting BAM Files @@@"
# Sort SAM files on chromosomal alignments
# http://www.htslib.org/doc/samtools.html
samtools sort -l 0 -o hisat2_sorted_221.bam -O bam hisat2_alignment_221.bam
samtools sort -l 0 -o hisat2_sorted_222.bam -O bam hisat2_alignment_222.bam

#echo "@@@ Retreiving Headers @@@"
# retreive the headers
#samtools view -H hisat2_sorted_217.bam

#echo "@@@ Running Cufflinks @@@"
# run Cufflinks
#cufflinks -o 217_clout_2 hisat2_sorted_217.bam
#cufflinks -o 219_clout_2 hisat2_sorted_219.bam

#echo "@@@ Running Cuffmerge  @@@"
#cuffmerge -g gencode.vM19.chr_patch_hapl_scaff.annotation.gtf assemblies.txt

#echo "@@@ Running Cuffdiff @@@"
#cuffdiff -o diff_out -u merged_asm/merged.gtf -L 217,219 -v  hisat2_sorted_217.bam hisat2_sorted_219.bam

### HTSeq (convert sorted files into count matrix for input into DESeq2
#htseq-count -f bam hisat2_sorted_217.bam hisat2_sorted_219.bam gencode.vM19.chr_patch_hapl_scaff.annotation.gtf

java -Xmx30G -jar /class/local/bin/QoRTs.jar QC --generatePlots --singleEnded --maxReadLength 50 ~/path/hisat2_sorted_221.bam ~/path/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf ~/path/qourts_output

R

install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source");
library(QoRTs)
res = read.qc.results.data("/users/adibe/Final-Project/QoRTs_Output/", decoder.files = "/users/adibe/Final-Project/Decoder/decoder.txt", calc.DESeq2 = TRUE)
makeMultiPlot.all(res, outfile.dir = "~/path/Summary_PDFs/", plot.device.name = "pdf")
quit()

java -Xmx30G -jar /class/local/bin/QoRTs.jar QC --generatePlots --singleEnded --maxReadLength 50 ~/path/hisat2_sorted_222.bam ~/path/gencode.vM19.chr_patch_hapl_scaff.annotation.gtf ~/path/qourts_output

R
##installing QoRTs package 
#install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos=NULL, type="source");
#library(QoRTs)

#Read in the QC data
res = read.qc.results.data("/users/adibe/Final-Project/QoRTs_Output/", decoder.files = "/users/adibe/Final-Project/Decoder/decoder.txt", calc.DESeq2 = TRUE)

#Create multiple pdfs for plots
makeMultiPlot.all(res, outfile.dir = "~/path/Summary_PDFs/", plot.device.name = "pdf")

#exit R with
quit()
