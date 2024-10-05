#!/bin/bash

######################### INPUT FILES #########################################
# Sample: Paired-end FASTQ files
# Reference Genome: hg38.fa
# Known Sites: dbSNP for hg38
# snpEff Database: GRCh38.99
###############################################################################

# Input FASTQ files (paired-end)
sample="SRR26688467_1.fastq SRR26688467_2.fastq"

# Reference genome
ref="hg38.fa"

# Paths to tools
picard="java -jar ~/NGS_course/tools/picard/picard.jar"
gatk="java -jar ~/NGS_course/tools/gatk/gatk-package-4.6.0.0-local.jar"
snpEff="java -jar ~/NGS_course/tools/snpEff/snpEff/snpEff.jar"
snpSift="java -jar ~/NGS_course/tools/snpEff/snpEff/SnpSift.jar"

# Known sites file (dbSNP)
knownsites="~/NGS_course/data/references/common_dbsnp.vcf.gz"

# snpEff database
snpeff_db="GRCh38.99"

###############################################################################
# Start of the shell script

# Set the date format
date_format="+%Y-%m-%d %H:%M:%S"

# Create necessary directories
mkdir -p ~/NGS_course/results/{qc_report,alignment,variants,annotations}

# Change to the raw data directory
cd ~/NGS_course/data/raw/

###############################################################################
# Step 1: Download data from SRA using fastq-dump (if not already downloaded)
###############################################################################

echo "$(date "$date_format") - Downloading FASTQ files"
fastq-dump --split-files SRR26688467

###############################################################################
# Step 2: Perform quality check with FastQC
###############################################################################

echo "$(date "$date_format") - Running FastQC"
fastqc -o ~/NGS_course/results/qc_report/ $sample

###############################################################################
# Step 3: Index the reference genome with BWA
###############################################################################

# Ensure the reference genome is uncompressed
echo "$(date "$date_format") - Preparing the reference genome"
cd ~/NGS_course/data/references/
if [ -f "${ref}.gz" ]; then
    gunzip "${ref}.gz"
fi

# Index the reference genome
bwa index $ref

###############################################################################
# Step 4: Align reads to the reference genome using BWA MEM
###############################################################################

echo "$(date "$date_format") - Aligning reads with BWA MEM"
cd ~/NGS_course/results/alignment/
bwa mem -M -t 6 -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tSM:sample_1' \
~/NGS_course/data/references/$ref \
~/NGS_course/data/raw/SRR26688467_1.fastq \
~/NGS_course/data/raw/SRR26688467_2.fastq \
> aligned_reads.sam

###############################################################################
# Step 5: Convert SAM to sorted BAM using Picard
###############################################################################

echo "$(date "$date_format") - Sorting SAM file and converting to BAM"
$picard SortSam \
    INPUT=aligned_reads.sam \
    OUTPUT=sorted_reads.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT

###############################################################################
# Step 6: Collect alignment and insert size metrics (optional)
###############################################################################

echo "$(date "$date_format") - Collecting alignment metrics"
$picard CollectAlignmentSummaryMetrics \
    R=~/NGS_course/data/references/$ref \
    I=sorted_reads.bam \
    O=alignment_metrics.txt

echo "$(date "$date_format") - Collecting insert size metrics"
$picard CollectInsertSizeMetrics \
    INPUT=sorted_reads.bam \
    OUTPUT=insert_metrics.txt \
    HISTOGRAM_FILE=insert_size_histogram.pdf

###############################################################################
# Step 7: Mark duplicates with Picard
###############################################################################

echo "$(date "$date_format") - Marking duplicates"
$picard MarkDuplicates \
    INPUT=sorted_reads.bam \
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

###############################################################################
# Step 8: Create sequence dictionary and index the reference genome
###############################################################################

echo "$(date "$date_format") - Creating sequence dictionary"
$picard CreateSequenceDictionary \
    R=~/NGS_course/data/references/$ref \
    O=~/NGS_course/data/references/${ref%.fa}.dict

echo "$(date "$date_format") - Indexing the reference genome with samtools"
samtools faidx ~/NGS_course/data/references/$ref

###############################################################################
# Step 9: Index the known sites VCF file
###############################################################################

echo "$(date "$date_format") - Indexing known sites VCF file"
cd ~/NGS_course/data/references/
if [ ! -f "${knownsites}.tbi" ]; then
    tabix -p vcf $knownsites
fi

###############################################################################
# Step 10: Base Quality Score Recalibration (BQSR) with GATK
###############################################################################

echo "$(date "$date_format") - Performing BaseRecalibrator"
cd ~/NGS_course/results/alignment/
$gatk BaseRecalibrator \
    -I dedup_reads.bam \
    -R ~/NGS_course/data/references/$ref \
    --known-sites $knownsites \
    -O recal_data.table

echo "$(date "$date_format") - Applying BQSR"
$gatk ApplyBQSR \
    -R ~/NGS_course/data/references/$ref \
    -I dedup_reads.bam \
    --bqsr-recal-file recal_data.table \
    -O recal_reads.bam

###############################################################################
# Step 11: Variant calling with GATK HaplotypeCaller
###############################################################################

echo "$(date "$date_format") - Calling variants with HaplotypeCaller"
$gatk HaplotypeCaller \
    -R ~/NGS_course/data/references/$ref \
    -I recal_reads.bam \
    -O ~/NGS_course/results/variants/raw_variants.vcf

###############################################################################
# Step 12: Select SNPs and INDELs from the VCF file
###############################################################################

echo "$(date "$date_format") - Selecting SNPs"
$gatk SelectVariants \
    -R ~/NGS_course/data/references/$ref \
    -V ~/NGS_course/results/variants/raw_variants.vcf \
    --select-type-to-include SNP \
    -O ~/NGS_course/results/variants/raw_snps.vcf

echo "$(date "$date_format") - Selecting INDELs"
$gatk SelectVariants \
    -R ~/NGS_course/data/references/$ref \
    -V ~/NGS_course/results/variants/raw_variants.vcf \
    --select-type-to-include INDEL \
    -O ~/NGS_course/results/variants/raw_indels.vcf

###############################################################################
# Step 13: Filter SNPs and INDELs based on quality metrics
###############################################################################

echo "$(date "$date_format") - Filtering SNPs"
$gatk VariantFiltration \
    -R ~/NGS_course/data/references/$ref \
    -V ~/NGS_course/results/variants/raw_snps.vcf \
    --filter-name "QD_filter" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS_filter" \
    --filter-expression "FS > 60.0" \
    --filter-name "MQ_filter" \
    --filter-expression "MQ < 40.0" \
    --filter-name "SOR_filter" \
    --filter-expression "SOR > 4.0" \
    -O ~/NGS_course/results/variants/filtered_snps.vcf

echo "$(date "$date_format") - Filtering INDELs"
$gatk VariantFiltration \
    -R ~/NGS_course/data/references/$ref \
    -V ~/NGS_course/results/variants/raw_indels.vcf \
    --filter-name "QD_filter" \
    --filter-expression "QD < 2.0" \
    --filter-name "FS_filter" \
    --filter-expression "FS > 200.0" \
    --filter-name "SOR_filter" \
    --filter-expression "SOR > 10.0" \
    -O ~/NGS_course/results/variants/filtered_indels.vcf

###############################################################################
# Step 14: Annotate variants using snpEff
###############################################################################

echo "$(date "$date_format") - Annotating SNPs with snpEff"
$snpEff -v $snpeff_db \
    ~/NGS_course/results/variants/filtered_snps.vcf \
    > ~/NGS_course/results/annotations/filtered_snps_ann.vcf

echo "$(date "$date_format") - Annotating INDELs with snpEff"
$snpEff -v $snpeff_db \
    ~/NGS_course/results/variants/filtered_indels.vcf \
    > ~/NGS_course/results/annotations/filtered_indels_ann.vcf

###############################################################################
# Step 15: Extract missense variants using SnpSift
###############################################################################

echo "$(date "$date_format") - Extracting missense variants"
$snpSift filter "(ANN[*].EFFECT = 'missense_variant')" \
    ~/NGS_course/results/annotations/filtered_snps_ann.vcf \
    > ~/NGS_course/results/annotations/missense_snps.vcf

###############################################################################
# End of the script
echo "$(date "$date_format") - Analysis complete"
