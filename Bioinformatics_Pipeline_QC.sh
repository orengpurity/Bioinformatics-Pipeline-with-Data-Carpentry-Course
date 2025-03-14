25 Feb 2025
## The pipeline
#Quality control - Assessing quality using FastQC
#Quality control - Trimming and/or filtering reads (if necessary)
#Align reads to reference genome
#Perform post-alignment clean-up
#Variant calling

## Transfer/download data to PC
cp -r /mnt/d/.backup/untrimmed_fastq
#downnloads zipped files of sequences

## Unzip files using gunzip
gunzip *.fastq.gz


### Quality Controls (view fastq files and clean/trim) ###

##  1)Quality Control Process
#assess the quality of the sequence reads contained in fastq files
#uses fastqc tool and produces fastq files 

# Details on the FASTQ format
#Line	Description
1	Always begins with ‘@’ and then information about the read
2	The actual DNA sequence
3	Always begins with a ‘+’ and sometimes the same info in line 1
4	Has a string of characters which represent the quality scores; must have same number of characters as line 2
#each sequence is made up of the 4 lines 
#Quality is interpreted as the probability of an incorrect base call
#Each character is assigned a quality score between 0 and 41 as shown in the chart below.

# Confirm presence of tools uses -h
fastqc -h

# Assessing quality using FastQC
#using a software program to assess read quality and filter out poor quality reads
#Has two steps:
#a) use FastQC program to visualize the quality of our reads (looks at quality collectively across all reads within a sample)
#b) use Trimmomatic program to filter out poor quality reads

# Fastqc plot interpretation
x-axis : bp positions and sequence length
y-axis : quality score
Green region : good quality scores (high quality score)
Yellow region : acceptable quality scores (1st to 3rd quartile range)
Red region : bad quality scores (median quality score)
Each bp position has a box-and-whisker plot that shows distribution of quality scores of all reads at a poition

# Step 1- Running FastQC
#FastQC can accept multiple file names as input, and on both zipped and unzipped files
#run FastQC on all of the FASTQ files in this directory
fastqc *.fastq*

# Step 2- Viewing Fastqc results
#For each input FASTQ file, FastQC creates a .zip file and a.html file. 
#The .html file is a stable webpage displaying the summary report for each of our samples
#The .html file is a stable webpage displaying the summary report for each of our samples
#move the results (.zip files and .html files into a results directory)

mkdir -p ~/dc_workshop/results/fastqc_untrimmed_reads
mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/

#view the .html files from PC a web browser (mozila forefox)
#If working on cloud, create a directory and transfer the .html files using scp command

## Decoding the FastQC Graphs Report (.html files)
#Per tile sequence quality; displays patterns in base quality along these tiles
#Per sequence quality scores; shows the most common quality scores
#Per base sequence content; proportion of each base position over all of the reads
#Per sequence GC content; average GC content in each of the reads
#Per base N content;  the percent of times that ‘aNy’ other nucleotides occurs at a position in all reads. Errors during sequencing
#Sequence Length Distribution; the distribution of sequence lengths of all reads in the file (its shorter in trimmed sequences)
#Sequence Duplication Levels; A distribution of duplicated sequences
#Overrepresented sequences; A list of sequences that occur more frequently than would be expected by chance
#Adapter Content; indicates where adapater sequences occur in the reads
#K-mer Content; showing any sequences which may show a positional bias within the reads

## Exploring FastQC text reports (.zip files)
#unzip/decompress the compressed files using unzip command
#unzip command unzips one file at a time; to unzip multiple files uses loop command
for filename in *.zip
> do
> unzip $filename
> done

## 2) Trimming/filter reads
#filter poor quality reads and trim poor quality bases of samples
#uses Trimmomatic tool and produces fastq files

#confirm  presnece of trimomatic
trimmoatic -h
#installing trimmomatic in ubuntu
Download miniconda
Install conda and activate
create conda environment
install trimmomatic by conda install -c bioconda trimmomatic

# Confirming trimommatic options to use
trimmomatic
# The options are either paired end (PE) or single end (SE) reads
#More trimmomatic options
step	meaning
ILLUMINACLIP	Perform adapter removal.
SLIDINGWINDOW	Perform sliding window trimming, cutting once the average quality within the window falls below a threshold.
LEADING	Cut bases off the start of a read, if below a threshold quality.
TRAILING	Cut bases off the end of a read, if below a threshold quality.
CROP	Cut the read to a specified length.
HEADCROP	Cut the specified number of bases from the start of the read.
MINLEN	Drop an entire read if it is below a specified length.
TOPHRED33	Convert quality scores to Phred-33.
TOPHRED64	Convert quality scores to Phred-64.
-threads 4	to use four computing threads to run (this will speed up our run)

# Step 1- Running Trimmomatic
#Trims zipped files
#running paired end samples (PE) for one sample 
trimmomatic PE SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz \
                SRR2589044_1.trim.fastq.gz SRR2589044_1un.trim.fastq.gz \
                SRR2589044_2.trim.fastq.gz SRR2589044_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
ls SRR2589044*
ls SRR2589044* -l -h
# using loop to run/trim multiple samples 
#zipping the unzipped files uses gzip command
gzip SRR2584863_1.fastq 
gzip  *.fastq
#looping
for infile in *_1.fastq.gz
> do
>   base=$(basename ${infile} _1.fastq.gz)
>   trimmomatic PE ${infile} ${base}_2.fastq.gz \
>                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
>                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
>                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
> done
#move the trimmed files to a new subdirectory
cd ~/dc_workshop/data/untrimmed_fastq
mkdir ../trimmed_fastq
mv *.trim* ../trimmed_fastq
cd ../trimmed_fastq
ls

## 3) Viewing trimmed sequences??? ##
fastqc ~/dc_workshop/data/trimmed_fastq/*.fastq*

### Next is : Alignment to a reference genome, cleaning of aligned sequence and variant calling
