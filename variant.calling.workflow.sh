12 March 2025

### Automating the variant calling workflow ###

## Downlading scripts ##
curl -O https://datacarpentry.org/wrangling-genomics/files/run_variant_calling.sh

## variant calling workflow ##
#The steps:
#Index the reference genome for use by bwa and samtools.
#Align reads to reference genome.
#Convert the format of the alignment to sorted BAM, with some intermediate steps.
#Calculate the read coverage of positions in the genome.
#Detect the single nucleotide variants (SNVs).
#Filter and report the SNVs in VCF (variant calling format)

## The Automation steps

#Preps.#
set -e for #the script to exit if an error occurs

#change directory
cd ~/dc_workshop/results

#find the reference genome by assigning the genome variable to the path to our reference genome
genome=~/dc_workshop/data/ref_genome/ecoli_rel606.fasta

# Step 1:index the reference genome for BWA
bwa index $genome
#create the directory structure to store the results
mkdir -p sam bam bcf vcf
#assign the name of the FASTQ file to a variable called fq1 while using loop (for all files)
for fq1 in ~/dc_workshop/data/trimmed_fastq_small/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"
#extract the base name of the file (excluding the path and .fastq extension) and assign it to a new variable called base
base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base" #echo is used within the command scripts to get an automated progress update#
#create variables to store the names of our output files (i.e both the base_1.fastq and base_2.fastq input files)
# add a different extension (e.g. .sam, .bam) for each file produced by our workflow

#input fastq files
    fq1=~/dc_workshop/data/trimmed_fastq_small/${base}_1.trim.sub.fastq
    fq2=~/dc_workshop/data/trimmed_fastq_small/${base}_2.trim.sub.fastq
 # output files
    sam=~/dc_workshop/results/sam/${base}.aligned.sam
    bam=~/dc_workshop/results/bam/${base}.aligned.bam
    sorted_bam=~/dc_workshop/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/dc_workshop/results/bcf/${base}_raw.bcf
    variants=~/dc_workshop/results/bcf/${base}_variants.vcf
    final_variants=~/dc_workshop/results/vcf/${base}_final_variants.vcf 

# Step 2: align the reads to the reference genome and output a .sam file
bwa mem $genome $fq1 $fq2 > $sam

# Step 3: convert the SAM file to BAM format
samtools view -S -b $sam > $bam

# Step 4:sort the BAM file
samtools sort -o $sorted_bam $bam

# Step 5: index the BAM file for display purposes
samtools index $sorted_bam

# Step 6: calculate the read coverage of positions in the genome
bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam 

# Step 7: call SNVs with bcftools
bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 

#Step 8: filter and report the SNVs in variant calling format (VCF)
vcfutils.pl varFilter $variants  > $final_variants

## Running the script
bash variant.calling.workflow.sh

