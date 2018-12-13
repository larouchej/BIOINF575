 #!/bin/bash

echo "@@@ Converting SRA to FASTQ @@@"
# SRA to FASTQ
cd tools/sratoolkit.2.9.2/bin
./fastq-dump -O ~/raw_fastq -p  SRR5590217
./fastq-dump -O ~/raw_fastq -p  SRR5590219
./fastq-dump -O ~/raw_fastq -p  SRR5590221
./fastq-dump -O ~/raw_fastq -p  SRR5590222

# Quality Control (using fastqc)
cd raw_fastq
fastqc SRR5590217.fastq SRR5590219.fastq SRR5590221.fastq SRR5590222.fastq

echo "@@@ Trimmomatic @@@"
# Trimmomatic (default parameters)
trim = tools/Trimmomatic-0.38/trimmomatic-0.38.jar
adapters = tools/Trimmomatic-0.38/adapters/TruSeq3-SE
java -jar $trim SE -phred33 raw_fastq/SRR5590221.fastq 02-trimmed_output/221.trimmed.fastq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar $trim SE -phred33 raw_fastq/SRR5590222.fastq 02-trimmed_output/222.trimmed.fastq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar $trim SE -phred33 raw_fastq/SRR5590217.fastq 02-trimmed_output/217.trimmed.fastq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
java -jar $trim SE -phred33 raw_fastq/SRR5590219.fastq 02-trimmed_output/219.trimmed.fastq ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Quality Control (again using fastqc)
fastqc 217.trimmed.fastq 219.trimmed.fastq 221.trimmed.fastq 222.trimmed.fastq

echo "@@@ Generating Kallisto index @@@"
kallisto index -i references/transcripts.idx ./Mus_musculus.GRCm38.cdna.all.fa.gz 

echo "@@@ Quantifying transcript abundances: SRR5590217 @@@"
 kallisto quant -i references/transcripts.idx \
 -g references/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz \
 -o 03-kallisto_output/217_Aged_Syn \
 -c references/chrom.sizes \z
 --genomebam --single -l 200 -s 20 -b 100 \
 02-trimmed_output/SRR5590217.trimmed.fastq.gz

echo "@@@ Quantifying transcript abundances: SRR5590219 @@@"
  kallisto quant -i references/transcripts.idx \
 -g references/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz \
 -o 03-kallisto_output/219_Adult_Syn \
 -c references/chrom.sizes \z
 --genomebam --single -l 200 -s 20 -b 100 \
 02-trimmed_output/SRR5590219.trimmed.fastq.gz

echo "@@@ Quantifying transcript abundances: SRR5590221 @@@"
  kallisto quant -i references/transcripts.idx \
 -g references/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz \
 -o 03-kallisto_output/221_Adult_Ctrl \
 -c references/chrom.sizes \z
 --genomebam --single -l 200 -s 20 -b 100 \
 02-trimmed_output/SRR5590221.trimmed.fastq.gz

echo "@@@ Quantifying transcript abundances: SRR5590222 @@@"
  kallisto quant -i references/transcripts.idx \
 -g references/gencode.vM6.chr_patch_hapl_scaff.annotation.gtf.gz \
 -o 03-kallisto_output/222_Aged_Ctrl \
 -c references/chrom.sizes \z
 --genomebam --single -l 200 -s 20 -b 100 \
 02-trimmed_output/SRR5590222.trimmed.fastq.gz