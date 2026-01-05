#!/bin/bash

# Script to filter and annotate variants
# This script is for demonstration purposes only


# directories
ref="/path/to/supporting_files/hg38/hg38.fa"
results="/path/to/results/SRR062634" 

if false
then

# -------------------
# Filter Variants - GATK4
# -------------------

# Filter SNPs
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/SRR062634_raw_snps.vcf \
	-O ${results}/SRR062634_filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



# Filter INDELS
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/SRR062634_raw_indels.vcf \
	-O ${results}/SRR062634_filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



# Select Variants that PASS filters
gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SRR062634_filtered_snps.vcf \
        -O ${results}/SRR062634_analysis-ready-snps.vcf


gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/SRR062634_filtered_indels.vcf \
	-O ${results}/SRR062634_analysis-ready-indels.vcf


# to exclude variants that failed genotype filters (cd results)
cat ${results}/SRR062634_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > ${results}/SRR062634_analysis-ready-snps-filteredGT.vcf
cat ${results}/SRR062634_analysis-ready-indels.vcf| grep -v -E "DP_filter|GQ_filter" > ${results}/SRR062634_analysis-ready-indels-filteredGT.vcf



# -------------------
# Annotate Variants - GATK4 Funcotator
# -------------------

# Annotate using Funcotator
gatk Funcotator \
	--variant ${results}/SRR062634_analysis-ready-snps-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /path/to/dataSources/funcotator_dataSources.v1.8.hg38.20230908g\
	--output ${results}/SRR062634_analysis-ready-snps-filteredGT-functotated.vcf \
	--output-file-format VCF

gatk Funcotator \
	--variant ${results}/SRR062634_analysis-ready-indels-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
        --data-sources-path /path/to/dataSources/funcotator_dataSources.v1.8.hg38.20230908g\
	--output ${results}/SRR062634_analysis-ready-indels-filteredGT-functotated.vcf \
	--output-file-format VCF


fi

# Extract fields from a VCF file to a tab-delimited table

gatk VariantsToTable \
	-V ${results}/SRR062634_analysis-ready-snps-filteredGT-functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/SRR062634_output_snps.table

# Extract information on a single gene, i.e. "NBPF1" in a clear excel-friendly table
cat ${results}/SRR062634_analysis-ready-snps-filteredGT-functotated.vcf|grep "Funcotation fields are: " |sed 's/|/\t/g'> output_NBPF1_curated_snps.txt
cat ${results}/SRR062634_output_snps.table|cut -f 5 |grep "NBPF1" |sed 's/|/\t/g'>> output_NBPF1_curated_snps.txt





