
26 Feb 2025

### Bioinformatics pipeline


### Alignment to a reference genome ###

#Tool used is Burrows Wheeler Aligner (BWA) and produces  SAM/BAM files
#BWA is a software package for mapping low-divergent sequences against a large reference genome
#Performs read alignment or mapping to determine where in the genome our reads originated from

## The alignment process consists of two steps:
# 1. Indexing the reference genome
# 2. Aligning the reads to the reference genome

#  Preps. #
#Download the reference genome
cd ~/dc_workshop
mkdir -p data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
#Unzip the reference genome
gunzip data/ref_genome/ecoli_rel606.fasta.gz
#Check real genaome's name
head data/ref_genome/ecoli_rel606.fasta

## Step 1 - Index the reference genome
#Indexing is used by BWA to quickly find potential alignment sites for query sequences in a genome
bwa index data/ref_genome/ecoli_rel606.fasta

## Step 2 - Align reads to reference genome
#It's the act of choosing an appropriate reference genome to map our reads against and then deciding on an aligner
#Uses the BWA-MEM algorithm (bwa mem ref_genome.fa, input_file_R1.fastq, or input_file_R2.fastq)
#Alignment is done using the trimmed sequneces from QC process

## Alingning a sample instance (SRR2584866)
#Move the trimed sequences to the current directory of reference genome (but should be in different sub directories)
bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam

## Alingnment output (SAM file)
#It's a tab-delimited text file that contains information for each individual read and its alignment to the genome

## Conversion of SAM to BAM
#BAM (reduced size and allows indexing) is the compressed binary version of SAM
#uses the samtools program with the view command and tell this command that the input is in SAM format (-S) and to output BAM format (-b)
samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam

## Sorting BAM file by coordinates
#Next we sort the BAM file using the sort command from samtools. -o tells the command where to write the output.
samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 

## Using samtools to learn more about this bam file
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam



### Next is Variant Calling Process ###


