10 March 2025

### Bioinformatics Pipeline

### Variant Calling ###

#A variant call is a conclusion that there is a nucleotide difference,
# vs. some reference at a given position in an individual genome or transcriptome, often referred to as a Single Nucleotide Variant (SNV)
#The call is accompanied by an estimate of variant frequency and some measure of confidence
#Tools used  are bcftools and produces a vcf file


## Preps. 

#Step 1: Calculate the read coverage of positions in the genome
#counting read coverage with bcftools
#uses the command mpileup (To generate a file with coverage information for every base)
#uses flags;
# -O b tells bcftools to generate a bcf format output file
# -o specifies where to write the output file
# -f flags the path to the reference genome
bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam

## Step 2: Detect the single nucleotide variants (SNVs)
#Identifing SNVs using bcftools call command
#Uses flags;
# --ploidy specify ploidy (which is one for the haploid E. coli)
# -m allows for multiallelic and rare-variant calling
# -v tells the program to output variant sites only (not every site in the genome)
# -o specifies where to write the output file
bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf

## Variant calling process ##

## Step 3: Filter and report the SNV variants in variant calling format (VCF)
#Uses the vcfutils.pl command to filter the SNVs for the final output in VCF format
vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf  > results/vcf/SRR2584866_final_variants.vcf

##Exploring the VCF format
less -S results/vcf/SRR2584866_final_variants.vcf
#VCF file infomation
# a) The header;
#the time and date the file was created
#the version of bcftools that was used
#the command line parameters used
# b) Additional informationon on each of the variations observed. They include:

# column	info
#CHROM	contig location where the variation occurs
#POS	position within the contig where the variation occurs
#ID	a . until we add annotation information
#REF	reference genotype (forward strand)
#ALT	sample genotype (forward strand)
#QUAL	Phred-scaled probability that the observed variant exists at this site (higher is better)
#FILTER	a . if no quality filters have been applied, PASS if a filter is passed, or the name of the filters this variant failed
#FORMAT	lists in order the metrics presented in the final column
#results	lists the values associated with those metrics in order

# c) Assessing how many variants are in the vcf file
grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l


## Assessing the alignment (visualization) ## (Optional)
#Visualizing the data in a genome browser to detect abnormalities and problems for further analysis
#Can be done using two visualization toola: 
# a)tview - a light-weight command-line based
# b)Institute’s Integrative Genomics Viewer (IGV) - requires software installation and transfer of files 

## Preps. 
# Step 1.  indexing the BAM file using samtools command
samtools index results/bam/SRR2584866.aligned.sorted.bam

## Visualizations

# Step 2(a). Viewing with tview
#Uses samtools command
#Uses tview, giving it the sorted bam file and the reference file
samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta

# Exploring the view
#The first line of output shows the genome coordinates in our reference genome
#The second line shows the reference genome sequence
#The third line shows the consensus sequence determined from the sequence reads
#The dot (.) indicates a match to the reference sequence

# Navigating to a specific position
#type g - A dialogue box will appear
#In this box, type the name of the “chromosome” followed by a colon and the position of the variant you would like to view
#For instance; for this sample, type CP000819.1:50 to view the 50th base

# 2(b). Viewing with IGV
#IGV is a stand-alone browser is installed locally and provide fast access
#Uses locally stored files (can use scp command to copy to local machine)

## Preps.
#Download IGV from the Broad Institute’s software page
#Open IGV
#Load our reference genome file (ecoli_rel606.fasta) into IGV using the “Load Genomes from File…” option under the “Genomes” pull-down menu.
#Load our BAM file (SRR2584866.aligned.sorted.bam) using the “Load from File…” option under the “File” pull-down menu.
#Do the same with our VCF file (SRR2584866_final_variants.vcf)

## Explore the view


### Automating a Variant Calling Workflow ###

# Uses the command touch to create a new file to write the shell script
# The command touch allows the creation of a new file without opening that file (unlike nano)

# Step 1. Create the file using touch
touch (file name)
ls

# Step 2. Open the file in nano and start building the script
nano (file name)
