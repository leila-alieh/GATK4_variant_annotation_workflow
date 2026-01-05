#!/bin/bash

# Script to call germline variants in a human WGS paired end reads 2 X 100bp: SRR062634.
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# With minor adaptation from https://github.com/kpatel427/YouTubeTutorials/blob/main/variant_calling.sh
# This script is for demonstration purposes only


# download data original
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

# download data: to avoid incomplete download and corrupted file use sra toolkit
prefetch SRR062634
fastq-dump --split-files path/to/sample_sra/SRR062634/SRR062634.sra -O /path/to/reads 

# compress reads
gzip /path/to/reads/SRR062634_1.fastq
gzip /path/to/reads/SRR062634_2.fastq

echo "Run Prep files..."

################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################

# download reference files
wget -P /path/to/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /path/to/supporting_files/hg38/hg38.fa.gz

# index ref - .fai file before running haplotype caller
samtools faidx /path/to/supporting_files/hg38/hg38.fa


# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=/path/to/supporting_files/hg38/hg38.fa O=/path/to/supporting_files/hg38/hg38.dict


# download known sites files for BQSR from GATK resource bundle
wget -P /path/to/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /path/to/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


###################################################### VARIANT CALLING STEPS ####################################################################


# directories
ref="/path/to/supporting_files/hg38/hg38.fa"
known_sites="/path/to/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/path/to/VC/aligned_reads"
reads="/path/to/VC/reads"
results="/path/to/VC/results"
data="/path/to/VC/data"


# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming required, quality looks okay.


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam




# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf






