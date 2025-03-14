12 March 2025

### Automating Scripts Workflow ####

## Analyzing quality with FastQC ##

set -e #(ensure that our script will exit if an error occurs)

cd ~/Data_Carpentry/scripts #(move me into the scripts/ directory when I run my script)

echo "Running FastQC ..."
fastqc *.fastq*

mkdir -p ~/dc_workshop/results/fastqc_untrimmed_reads

echo "Saving FastQC results..."
mv *.zip ~/dc_workshop/results/fastqc_untrimmed_reads/
mv *.html ~/dc_workshop/results/fastqc_untrimmed_reads/
cd ~/dc_workshop/results/fastqc_untrimmed_reads/

echo "Unzipping..."
for filename in *.zip
    do
    unzip $filename
    done

echo "Saving summary..."
cat */summary.txt > ~/dc_workshop/docs/fastqc_summaries.txt

#Running the script
bash read_qc.sh

#Running QC for other samples using the sampe script#

#Select yes (after the first run is done) promt to replace the first sample with the next one (and so forth)
replace SRR2584866_fastqc/Icons/fastqc_icon.png? [y]es, [n]o, [A]ll, [N]one, [r]ename:
A













