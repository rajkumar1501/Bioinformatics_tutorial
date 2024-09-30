# Next-Generation Sequencing (NGS) Analysis Tutorial on Linux

Welcome to this comprehensive tutorial on performing Next-Generation Sequencing (NGS) analysis using a newly built Linux setup. This guide is designed for students and researchers new to NGS and bioinformatics. It provides step-by-step instructions on setting up the required applications, organizing your workspace, and conducting NGS analysis using real data.

We will use the Autism Multiplex data from Indian Families (SRR26688467) as a case study, guiding you through the entire process from system setup to variant annotation.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Directory Structure Setup](#directory-structure-setup)
3. [System Setup](#system-setup)
4. [Installing Required Applications](#installing-required-applications)
5. [Downloading Reference Genomes and Databases](#downloading-reference-genomes-and-databases)
6. [Data Acquisition](#data-acquisition)
7. [Quality Control with FastQC](#quality-control-with-fastqc)
8. [Sequence Alignment with BWA](#sequence-alignment-with-bwa)
9. [Sorting and Processing Alignments](#sorting-and-processing-alignments)
10. [Variant Calling and Annotation](#variant-calling-and-annotation)
    - [Installing and Running snpEff](#installing-and-running-snpeff)
11. [Summary](#summary)
12. [Additional Resources](#additional-resources)

---

## Introduction

### Study Overview

- **Study:** Whole Exome Sequencing on a multiplex family of Indian origin
- **SRX ID:** SRX22388233
- **SRA Run ID:** SRR26688467
- **Instrument:** Illumina HiSeq X
- **Strategy:** Whole Exome Sequencing (WXS)
- **Organism:** *Homo sapiens*
- **Sample:** 01S2

**Data Summary:**

- **Run:** 1 run
- **Number of Spots:** 23.5 Million
- **Number of Bases:** 6.7 Gigabases
- **Download Size:** 2.8 Gigabytes
- **Read Length:** 2x150 bp
- **Sequencing Depth:** 80-100X on target

### Experimental Design

Whole-Exome Sequencing (WXS) libraries were prepared using the SureSelectXT Human All Exon (V5) kit. The workflow included DNA shearing, end repair, adenylation, adapter ligation, and hybridization with exome-specific biotinylated capture probes. The enriched libraries were sequenced on the Illumina HiSeq X platform to generate paired-end reads.

---

## Directory Structure Setup

### Importance of an Organized Directory Structure

An organized directory structure is crucial for efficient data management, reproducibility, and collaboration. It helps in keeping raw data, processed data, tools, scripts, and results systematically arranged, making the analysis workflow more streamlined.

### Recommended Directory Layout

We recommend the following directory structure for your NGS analysis project:

```
~/NGS_course/
├── data/
│   ├── raw/
│   ├── processed/
│   └── references/
├── tools/
│   ├── bwa/
│   ├── samtools/
│   ├── fastqc/
│   ├── picard/
│   ├── snpEff/
│   ├── gatk/
│   └── sra-toolkit/
├── results/
│   ├── qc_report/
│   ├── alignment/
│   ├── variants/
│   └── annotations/
├── scripts/
└── logs/
```

**Directory Breakdown:**

- **data/**: Contains all data-related files.
  - **raw/**: Stores raw sequencing data (FASTQ files).
  - **processed/**: Holds processed data such as aligned reads (SAM/BAM files) and variant calls (VCF files).
  - **references/**: Contains reference genomes, annotation files, and databases.

- **tools/**: Houses all bioinformatics tools and their related files.
  - Each tool (e.g., BWA, Samtools) has its own subdirectory for organization.

- **results/**: Stores output files from various analysis steps.
  - **qc_report/**: FastQC reports.
  - **alignment/**: Aligned SAM/BAM files.
  - **variants/**: Variant call files (VCF).
  - **annotations/**: Annotated variant files.

- **scripts/**: Contains custom scripts for automating tasks.

- **logs/**: Stores log files generated during analysis.

### Creating the Directory Structure

**Command:**

```bash
# Navigate to home directory
cd ~

# Create main project directory with subdirectories
mkdir -p NGS_course/{data/{raw,processed,references},tools/{bwa,samtools,fastqc,picard,snpEff,gatk,sra-toolkit},results/{qc_report,alignment,variants,annotations},scripts,logs}

# Verify the directory structure
tree NGS_course
```

**Rationale:**

- **`mkdir -p`:** Creates directories and parent directories as needed.
- **`tree NGS_course`:** Displays the directory structure to verify correctness (install `tree` if necessary).

---

## System Setup

### Update the System

It's essential to keep your system up-to-date to ensure compatibility and security.

**Commands:**

```bash
sudo apt-get update
sudo apt-get upgrade -y
```

**Rationale:**

- **`sudo apt-get update`:** Updates the list of available packages.
- **`sudo apt-get upgrade -y`:** Installs the newest versions of all packages.

### Verify Current Directory

Before proceeding, ensure you're in the correct working directory.

**Commands:**

```bash
pwd
ls
```

**Rationale:**

- **`pwd`:** Displays the current working directory.
- **`ls`:** Lists the contents of the directory.

### Navigate to the Project Directory

**Command:**

```bash
cd ~/NGS_course/
```

**Rationale:**

- Moves you into the main project directory where all subsequent work will be conducted.

---

## Installing Required Applications

We will install the necessary bioinformatics tools within the `tools/` directory to keep them organized and separate from system-wide installations.

### General Installation Note

Some tools will be installed using `apt-get` (system-wide), while others will be downloaded and placed in specific directories.

### 1. BWA (Burrows-Wheeler Aligner)

**Purpose:** Align sequencing reads to a reference genome.

**Commands:**

```bash
# Navigate to tools directory
cd ~/NGS_course/tools/bwa/

# Install BWA using apt-get
sudo apt-get install bwa -y

# Verify installation
bwa
```

**Rationale:**

- **`sudo apt-get install bwa -y`:** Installs BWA.
- **`bwa`:** Checks that BWA is accessible from the command line.

### 2. Samtools

**Purpose:** Manipulate SAM/BAM files.

**Commands:**

```bash
# Navigate to tools directory
cd ~/NGS_course/tools/samtools/

# Install Samtools using apt-get
sudo apt-get install samtools -y

# Verify installation
samtools --version
```

**Rationale:**

- Ensures that Samtools is installed and its version is displayed.

### 3. FastQC

**Purpose:** Perform quality control checks on raw sequencing data.

**Commands:**

```bash
# Navigate to tools directory
cd ~/NGS_course/tools/fastqc/

# Install FastQC using apt-get
sudo apt-get install fastqc -y

# Verify installation
fastqc --version
```

**Rationale:**

- Confirms FastQC is installed and ready for use.

### 4. Tabix

**Purpose:** Index genomic data in tab-delimited files.

**Commands:**

```bash
# Navigate to tools directory
cd ~/NGS_course/tools/tabix/

# Install Tabix using apt-get
sudo apt-get install tabix -y

# Verify installation
tabix --version
```

**Rationale:**

- Ensures Tabix is installed correctly.

### 5. Picard Tools

**Purpose:** Manipulate high-throughput sequencing (HTS) data and formats.

**Commands:**

```bash
# Navigate to Picard tools directory
cd ~/NGS_course/tools/picard/

# Download Picard JAR file
wget https://github.com/broadinstitute/picard/releases/download/2.18.1/picard.jar

# Verify download
ls picard.jar
```

**Rationale:**

- Downloads Picard directly into the tools directory.
- Verifies that `picard.jar` is present.

### 6. snpEff

**Purpose:** Annotate and predict effects of genetic variants.

#### Installing snpEff

**Commands:**

```bash
# Navigate to tools directory
cd ~/NGS_course/tools/snpEff/

# Download snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip snpEff
unzip snpEff_latest_core.zip

# Navigate to snpEff directory
cd snpEff/

# Verify installation
ls snpEff.jar
```

**Rationale:**

- Downloads and unzips snpEff in the appropriate directory.
- Ensures that the `snpEff.jar` file is present.

#### Downloading snpEff Databases

Since we are using the hg38 reference genome, we need to download the corresponding snpEff database `GRCh38.99`.

**Commands:**

```bash
# Download the database for GRCh38.99
java -jar snpEff.jar download GRCh38.99

# Verify available databases
java -jar snpEff.jar databases | grep -i GRCh38
```

**Rationale:**

- Downloads the database corresponding to the GRCh38 (hg38) reference genome.
- Verifies that the database is installed.

### 7. GATK (Genome Analysis Toolkit)

**Purpose:** Perform variant discovery in high-throughput sequencing data.

**Commands:**

```bash
# Navigate to GATK tools directory
cd ~/NGS_course/tools/gatk/

# Download GATK
wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip

# Unzip GATK
unzip gatk-4.6.0.0.zip

# Move GATK jar file to tools directory
mv gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar ./

# Verify installation
java -jar gatk-package-4.6.0.0-local.jar --help
```

**Rationale:**

- Downloads GATK and places the JAR file in the tools directory.
- Checks that GATK is accessible.

### 8. SRA Toolkit

**Purpose:** Download and convert SRA data to FASTQ format.

**Commands:**

```bash
# Navigate to SRA Toolkit tools directory
cd ~/NGS_course/tools/sra-toolkit/

# Download SRA Toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

# Extract the toolkit
tar -xzf sratoolkit.current-ubuntu64.tar.gz

# Add SRA Toolkit to PATH
echo 'export PATH=$PATH:~/NGS_course/tools/sra-toolkit/sratoolkit.*/bin' >> ~/.bashrc
source ~/.bashrc

# Verify installation
fastq-dump --version
```

**Rationale:**

- Downloads and extracts the SRA Toolkit.
- Adds the toolkit's `bin` directory to the PATH variable for easy access.
- Verifies that `fastq-dump` is available.

### 9. Additional Utilities

**Commands:**

```bash
sudo apt-get install libxml2-utils -y
```

**Rationale:**

- Installs utilities required by some bioinformatics tools (e.g., `xmllint`).

---

## Downloading Reference Genomes and Databases

### 1. Human Reference Genome (hg38)

**Purpose:** Provide a reference genome for read alignment.

**Commands:**

```bash
# Navigate to references directory
cd ~/NGS_course/data/references/

# Download the hg38 reference genome
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Verify download
ls -lh hg38.fa.gz
```

**Rationale:**

- Downloads the compressed FASTA file for the hg38 human reference genome.
- Ensures the file is present and checks its size.

### 2. dbSNP Database

**Purpose:** Annotate known variants during variant calling and annotation.

**Commands:**

```bash
# Navigate to references directory
cd ~/NGS_course/data/references/

# Download dbSNP VCF file
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/common_all_20180418.vcf.gz

# Rename for clarity
mv common_all_20180418.vcf.gz common_dbsnp.vcf.gz

# Verify download
ls -lh common_dbsnp.vcf.gz
```

**Rationale:**

- Downloads the dbSNP database in VCF format.
- Renames the file for easier identification.
- Confirms the presence of the file.

---

## Data Acquisition

### Download FASTQ Files Using SRA Toolkit

**Purpose:** Obtain the raw sequencing data for analysis.

**Commands:**

```bash
# Navigate to raw data directory
cd ~/NGS_course/data/raw/

# Download FASTQ files
fastq-dump --split-files SRR26688467

# List the downloaded files
ls -lh SRR26688467_*
```

**Rationale:**

- **`fastq-dump --split-files`:** Downloads and splits paired-end reads into separate files.
- Ensures that the FASTQ files (`SRR26688467_1.fastq` and `SRR26688467_2.fastq`) are present and checks their sizes.

---

## Quality Control with FastQC

### Purpose

Assess the quality of the raw sequencing data to identify any issues before proceeding with analysis.

### Steps

1. **Create a QC Report Directory**

   **Command:**

   ```bash
   mkdir -p ~/NGS_course/results/qc_report/
   ```

   **Rationale:**

   - Organizes QC reports in a dedicated directory.

2. **Run FastQC**

   **Command:**

   ```bash
   fastqc -o ~/NGS_course/results/qc_report/ ~/NGS_course/data/raw/SRR26688467_1.fastq ~/NGS_course/data/raw/SRR26688467_2.fastq
   ```

   **Rationale:**

   - **`-o`:** Specifies the output directory.
   - Processes both forward and reverse reads.

3. **View QC Reports**

   **Commands:**

   ```bash
   cd ~/NGS_course/results/qc_report/
   ls
   ```

   **Rationale:**

   - Lists the generated QC reports.
   - Reports can be viewed using a web browser to assess data quality metrics like base quality scores, GC content, and sequence duplication levels.

---

## Sequence Alignment with BWA

### Purpose

Align sequencing reads to the reference genome to identify where each read originated from.

### Steps

1. **Index the Reference Genome**

   **Commands:**

   ```bash
   cd ~/NGS_course/data/references/

   # Decompress the reference genome if necessary
   gunzip hg38.fa.gz

   # Index the reference genome
   bwa index hg38.fa

   # Verify index files
   ls hg38.fa.*
   ```

   **Rationale:**

   - Decompresses the reference genome if it's compressed (BWA requires an uncompressed FASTA file).
   - Prepares the reference genome for alignment.
   - BWA creates index files required for efficient alignment.

2. **Perform Alignment**

   **Create Alignment Directory**

   ```bash
   mkdir -p ~/NGS_course/results/alignment/
   ```

   **Run BWA MEM**

   ```bash
   bwa mem -M -t 6 -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tSM:sample_1' \
   ~/NGS_course/data/references/hg38.fa \
   ~/NGS_course/data/raw/SRR26688467_1.fastq \
   ~/NGS_course/data/raw/SRR26688467_2.fastq \
   > ~/NGS_course/results/alignment/aligned_reads.sam
   ```

   **Rationale:**

   - **`-M`:** Marks shorter split hits as secondary, important for downstream compatibility with GATK.
   - **`-t 6`:** Uses 6 threads for faster processing.
   - **`-R`:** Adds read group information, crucial for downstream analyses like variant calling.
     - **ID:** Read group identifier.
     - **LB:** Library identifier.
     - **PL:** Platform/technology used (e.g., ILLUMINA).
     - **SM:** Sample name.
   - Aligns paired-end reads to the reference genome.
   - Outputs the alignment in SAM format.

---

## Sorting and Processing Alignments

### Purpose

Convert SAM files to BAM format, sort them, and prepare for variant calling.

### Steps

1. **Sort the SAM File**

   **Commands:**

   ```bash
   cd ~/NGS_course/results/alignment/

   # Convert SAM to sorted BAM using Picard
   java -jar ~/NGS_course/tools/picard/picard.jar SortSam \
   INPUT=aligned_reads.sam \
   OUTPUT=sorted_reads.bam \
   SORT_ORDER=coordinate \
   VALIDATION_STRINGENCY=SILENT
   ```

   **Rationale:**

   - Converts the large SAM file to a compressed BAM file.
   - Sorts the reads by coordinate, which is required for many downstream tools.
   - **`VALIDATION_STRINGENCY=SILENT`:** Suppresses validation errors (useful if your SAM file has minor formatting issues).

2. **Verify the Sorted BAM File**

   **Command:**

   ```bash
   ls -lh sorted_reads.bam
   ```

   **Rationale:**

   - Checks the size and presence of the sorted BAM file.

3. **Index the Sorted BAM File**

   **Command:**

   ```bash
   samtools index sorted_reads.bam
   ```

   **Rationale:**

   - Creates an index file for the BAM file, allowing for fast random access during variant calling.

---

## Variant Calling and Annotation

### Purpose

Identify genetic variants (e.g., SNPs, indels) from the aligned reads and annotate them for biological significance.

### Steps

1. **Optional: Mark Duplicates**

   **Commands:**

   ```bash
   # Mark duplicates using Picard
   java -jar ~/NGS_course/tools/picard/picard.jar MarkDuplicates \
   INPUT=sorted_reads.bam \
   OUTPUT=dedup_reads.bam \
   METRICS_FILE=dedup_metrics.txt \
   CREATE_INDEX=true \
   VALIDATION_STRINGENCY=SILENT

   # Verify the deduplicated BAM file
   ls -lh dedup_reads.bam
   ```

   **Rationale:**

   - Identifies and marks duplicate reads that may arise from PCR amplification.
   - Reduces false-positive variant calls.
   - **`CREATE_INDEX=true`:** Automatically indexes the deduplicated BAM file.
   - **`METRICS_FILE`:** Outputs duplication metrics.

2. **Variant Calling with GATK**

   **Commands:**

   ```bash
   # Navigate to alignment directory
   cd ~/NGS_course/results/alignment/

   # Call variants using GATK HaplotypeCaller
   java -jar ~/NGS_course/tools/gatk/gatk-package-4.6.0.0-local.jar HaplotypeCaller \
   -R ~/NGS_course/data/references/hg38.fa \
   -I dedup_reads.bam \
   -O ~/NGS_course/results/variants/raw_variants.vcf.gz \
   -ERC GVCF \
   --native-pair-hmm-threads 6
   ```

   **Rationale:**

   - **`-ERC GVCF`:** Outputs a GVCF file, useful for joint genotyping if you have multiple samples.
   - **`--native-pair-hmm-threads`:** Allocates threads to speed up computation.
   - Outputs a VCF file containing raw variant calls.

3. **Annotate Variants with snpEff**

   **Refer to the [Installing and Running snpEff](#installing-and-running-snpeff) section for detailed instructions.**

---

## Installing and Running snpEff

### Purpose

Annotate variants to predict their impact on genes and proteins using the correct reference genome database.

### Steps

1. **Navigate to snpEff Directory**

   **Command:**

   ```bash
   cd ~/NGS_course/tools/snpEff/snpEff/
   ```

   **Rationale:**

   - Ensures you're in the correct directory to run snpEff.

2. **Run snpEff**

   **Commands:**

   ```bash
   java -Xmx8g -jar snpEff.jar GRCh38.99 \
   ~/NGS_course/results/variants/raw_variants.vcf.gz \
   > ~/NGS_course/results/annotations/annotated_variants.vcf
   ```

   **Rationale:**

   - **`-Xmx8g`:** Allocates 8GB of memory.
   - **`GRCh38.99`:** Specifies the genome build matching hg38.
   - Outputs an annotated VCF file.

3. **Verbose Mode (Optional)**

   **Command:**

   ```bash
   java -Xmx8g -jar snpEff.jar -v GRCh38.99 \
   ~/NGS_course/results/variants/raw_variants.vcf.gz \
   > ~/NGS_course/results/annotations/annotated_variants.vcf
   ```

   **Rationale:**

   - **`-v`:** Provides detailed logging for troubleshooting.

4. **View Annotated Variants**

   **Command:**

   ```bash
   head ~/NGS_course/results/annotations/annotated_variants.vcf
   ```

   **Rationale:**

   - Checks the annotated VCF file to confirm that annotations are present.

---

## Summary

### What We've Achieved

- **System and Directory Setup:** Established a well-organized directory structure and updated the Linux system.
- **Installed Essential Tools:** Installed BWA, Samtools, FastQC, Tabix, Picard, snpEff, GATK, and SRA Toolkit in a structured manner.
- **Downloaded References and Databases:** Acquired the hg38 reference genome and dbSNP database.
- **Data Acquisition:** Downloaded raw sequencing data using the SRA Toolkit.
- **Quality Control:** Assessed the quality of raw data with FastQC.
- **Read Alignment:** Aligned reads to the reference genome using BWA.
- **Data Processing:** Converted and sorted alignment files using Picard and Samtools.
- **Variant Calling:** Identified genetic variants using GATK HaplotypeCaller.
- **Variant Annotation:** Annotated variants with snpEff using the correct database for hg38 (GRCh38.99).

### Next Steps

- **Variant Filtering:** Apply filters to identify high-confidence variants.
- **Interpretation:** Analyze the biological significance of the variants.
- **Reporting:** Generate reports and visualizations for presentation.
- **Automation:** Develop scripts to automate the pipeline for future analyses.

---

## Additional Resources

- **BWA Documentation:** [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)
- **Samtools Documentation:** [http://www.htslib.org/doc/samtools.html](http://www.htslib.org/doc/samtools.html)
- **FastQC Documentation:** [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **Picard Tools:** [https://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)
- **snpEff Documentation:** [http://snpeff.sourceforge.net/](http://snpeff.sourceforge.net/)
- **GATK Documentation:** [https://gatk.broadinstitute.org/hc/en-us](https://gatk.broadinstitute.org/hc/en-us)
- **SRA Toolkit Documentation:** [https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)

---

**Note:** Always ensure that you are using compatible versions of tools and reference genomes to avoid inconsistencies. For instance, when using the hg38 reference genome, make sure to download the corresponding snpEff database (`GRCh38.99`) to match the genome build. Regularly update your tools and databases to utilize the latest features and improvements.