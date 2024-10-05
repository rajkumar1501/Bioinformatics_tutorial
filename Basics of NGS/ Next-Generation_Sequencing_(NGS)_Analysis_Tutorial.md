# Comprehensive NGS Data Analysis Tutorial for Beginners

Welcome to this detailed tutorial on Next-Generation Sequencing (NGS) data analysis using a Linux system. This guide is tailored for beginners and covers each step extensively, explaining the purpose, command parameters, alternative tools, and relevant information to help you understand and perform NGS data analysis effectively.

We will use a publicly available dataset (SRA accession number: **SRR26688467**) from a study on autism in Indian families. This dataset represents whole exome sequencing data from an individual with autism, providing a real-world example for our analysis.

---

## Introduction

### Study Overview

- **Study:** Whole Exome Sequencing of an individual with autism from an Indian multiplex family.
- **SRA Run ID:** SRR26688467
- **Instrument:** Illumina HiSeq X
- **Strategy:** Whole Exome Sequencing (WXS)
- **Organism:** *Homo sapiens*
- **Sample:** 01S2

**Data Summary:**

- **Number of Reads:** Approximately 23.5 million paired-end reads
- **Read Length:** 2x150 bp
- **Sequencing Depth:** 80-100X on target regions

### Objective

The goal of this tutorial is to guide you through the entire NGS data analysis pipeline, from raw data acquisition to variant annotation, providing detailed explanations suitable for beginners.

---

## Prerequisites

### Hardware Requirements

- A Linux system with at least 8 GB of RAM (16 GB or more recommended)
- Sufficient storage space (at least 50 GB) for data and intermediate files

### Software Requirements

Ensure that the following software is installed on your system:

- **Java (JDK 1.8 or higher)**
- **BWA (Burrows-Wheeler Aligner)**
- **Samtools**
- **FastQC**
- **Picard Tools**
- **GATK (Genome Analysis Toolkit) version 4.x**
- **snpEff**
- **SRA Toolkit**
- **Tabix**
- **SnpSift** (part of snpEff)

### Installation Resources
Before beginning the NGS analysis, ensure that your Linux system is up-to-date and properly configured.

1. **Update the System:**

   ```bash
   sudo apt-get update
   sudo apt-get upgrade -y
   ```

2. **Verify the Current Directory:**

   ```bash
   pwd
   ls
   ```

3. **Create a Working Directory for the NGS Course:**

   ```bash
   mkdir ~/NGS_course
   cd ~/NGS_course/
   ```

---

## Installing Required Applications

Several bioinformatics tools are essential for NGS analysis. Below are the installation steps for each tool.

### 1. **BWA (Burrows-Wheeler Aligner)**

BWA is used for aligning sequencing reads to a reference genome.

```bash
sudo apt-get install bwa -y
```

### 2. **Samtools**

Samtools is used for manipulating SAM/BAM files.

```bash
sudo apt-get install samtools -y
```

### 3. **FastQC**

FastQC is a quality control tool for high-throughput sequence data.

```bash
sudo apt-get install fastqc -y
```

### 4. **Tabix**

Tabix indexes genomic data in TAB-delimited files.

```bash
sudo apt-get install tabix -y
```

### 5. **Picard Tools**

Picard is a set of command-line tools for manipulating high-throughput sequencing (HTS) data and formats.

- **Download Picard:**

  ```bash
  wget https://github.com/broadinstitute/picard/releases/download/2.18.1/picard.jar
  ```

### 6. **snpEff**

snpEff is a tool for annotating and predicting the effects of genetic variants.

- **Download snpEff:**

  ```bash
  wget https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_GRCh37.75.zip
  ```

- **Unzip snpEff:**

  ```bash
  unzip snpEff_v5_0_GRCh37.75.zip
  cd snpEff_latest_core/snpEff/
  ```

- **Verify Available Databases:**

  ```bash
  java -jar snpEff.jar databases | grep -i homo_sapiens
  ```

### 7. **GATK (Genome Analysis Toolkit)**

GATK is used for variant discovery in high-throughput sequencing data.

- **Download GATK:**

  ```bash
  wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip
  unzip gatk-4.6.0.0.zip
  ```

- **Move GATK Jar to Working Directory:**

  ```bash
  mv gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar ~/NGS_course/
  ```

### 8. **sra-toolkit**

The SRA Toolkit is used for downloading and converting SRA data to FASTQ format.

- **Download SRA Toolkit:**

  ```bash
  wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
  tar -xzf sratoolkit.current-ubuntu64.tar.gz
  cd sratoolkit.3.1.1-ubuntu64/
  ```

- **Add SRA Toolkit to PATH:**

  ```bash
  echo "export PATH=$(pwd)/bin:\$PATH" >> ~/.bashrc
  source ~/.bashrc
  ```

### 9. **Additional Utilities**

Install additional necessary utilities:

```bash
sudo apt-get install libxml2-utils -y
```

---

## Downloading Reference Genomes and Databases

### 1. **Human Reference Genome (hg38)**

Download the human reference genome (GRCh38/hg38):

```bash
cd ~/NGS_course/
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```

- **Verify Download and Cleanup:**

  ```bash
  ls
  du -sh hg38.fa.gz
  ```

### 2. **dbSNP Database**

Download the dbSNP VCF file for variant annotation:

```bash
cd ~/NGS_course/snpEff_latest_core/snpEff/
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/common_all_20180418.vcf.gz
mv common_all_20180418.vcf.gz common_dbsnp.vcf.gz
```

---



## Setting Up the Working Environment

### Directory Structure

Creating a well-organized directory structure is crucial for managing your data and results effectively.

```bash
mkdir -p ~/NGS_tutorial/{data/{raw,processed,references},results/{qc_report,alignment,variants,annotations},tools,scripts,logs}
```

**Explanation:**

- **`~/NGS_tutorial/`**: Main project directory.
- **`data/`**: Contains raw and processed data.
  - **`raw/`**: Stores raw sequencing data (FASTQ files).
  - **`processed/`**: Stores processed data (BAM, VCF files).
  - **`references/`**: Contains reference genomes and databases.
- **`results/`**: Stores output from various analysis steps.
  - **`qc_report/`**: FastQC reports.
  - **`alignment/`**: Alignment files (SAM/BAM).
  - **`variants/`**: Variant files (VCF).
  - **`annotations/`**: Annotated variant files.
- **`tools/`**: Place for installing tools if needed.
- **`scripts/`**: Contains scripts for automation.
- **`logs/`**: Stores log files.

---

## Overview of the Analysis Pipeline

1. **Data Acquisition:** Download raw sequencing data from SRA.
2. **Quality Assessment:** Evaluate raw data quality using FastQC.
3. **Reference Genome Preparation:** Download and prepare the reference genome.
4. **Alignment:** Map reads to the reference genome using BWA MEM.
5. **Post-Alignment Processing:** Sort, convert, and mark duplicates.
6. **Base Quality Score Recalibration (BQSR):** Correct systematic errors in base quality scores.
7. **Variant Calling:** Identify genetic variants using GATK HaplotypeCaller.
8. **Variant Filtering:** Apply quality filters to variants.
9. **Variant Annotation:** Annotate variants with snpEff.
10. **Variant Extraction:** Extract specific variants of interest using SnpSift.

---

## Detailed Step-by-Step Tutorial

### Step 1: Data Acquisition

**Objective:** Obtain raw sequencing data for analysis.

**Command:**

```bash
cd ~/NGS_tutorial/data/raw/
fastq-dump --split-files SRR26688467
```

**Explanation:**

- **`cd ~/NGS_tutorial/data/raw/`**: Navigate to the directory where raw data will be stored.
- **`fastq-dump --split-files SRR26688467`**: Downloads the data and splits paired-end reads into two files:
  - **`SRR26688467_1.fastq`**: Forward reads.
  - **`SRR26688467_2.fastq`**: Reverse reads.

**Parameters:**

- **`--split-files`**: Splits paired-end reads into separate files.
- **`SRR26688467`**: SRA accession number for the sample.

**Why This Step is Important:**

- Acquiring raw data is essential to start the analysis.

**Alternative Tools:**

- **`prefetch`** and **`fasterq-dump`** from the SRA Toolkit can be used for faster downloads and conversions.

---

### Step 2: Quality Assessment with FastQC

**Objective:** Assess the quality of raw sequencing data.

**Command:**

```bash
mkdir -p ~/NGS_tutorial/results/qc_report/
fastqc -o ~/NGS_tutorial/results/qc_report/ SRR26688467_1.fastq SRR26688467_2.fastq
```

**Explanation:**

- **`fastqc`**: Runs quality checks on the FASTQ files.
- **`-o`**: Specifies the output directory for the reports.

**Parameters:**

- **`SRR26688467_1.fastq` and `SRR26688467_2.fastq`**: Input FASTQ files.
- **`-o ~/NGS_tutorial/results/qc_report/`**: Output directory for the FastQC reports.

**Why This Step is Important:**

- Identifies potential issues such as low-quality scores, adapter contamination, GC bias, and overrepresented sequences.
- Helps decide if data trimming or cleaning is required.

**Alternative Tools:**

- **FastP**: Performs both quality control and trimming.
- **MultiQC**: Aggregates results from multiple samples.

---

### Step 3: Preparing the Reference Genome

**Objective:** Download and prepare the reference genome for alignment and variant calling.

**Command:**

```bash
cd ~/NGS_tutorial/data/references/
# Download the hg38 reference genome
gunzip hg38.fa.gz
```

**Explanation:**

- **`wget`**: Downloads files from the internet.
- **`-O hg38.fa.gz`**: Saves the downloaded file as `hg38.fa.gz`.
- **`gunzip hg38.fa.gz`**: Decompresses the gzipped reference genome.

**Parameters:**

- **URL**: The link to the hg38 reference genome.
  - **Note**: Ensure you have the correct version matching your analysis.

**Why This Step is Important:**

- A reference genome is required for aligning reads and identifying variants.

**Alternative Tools:**

- **`rsync`**: Can be used for downloading files.
- **Reference Genomes from UCSC or Ensembl**: Alternative sources for genomes.

---

**Indexing the Reference Genome with BWA**

**Command:**

```bash
bwa index hg38.fa
```

**Explanation:**

- **`bwa index hg38.fa`**: Indexes the reference genome to prepare for alignment.

**Parameters:**

- **`hg38.fa`**: The reference genome FASTA file.

**Why This Step is Important:**

- Indexing allows BWA to perform efficient alignments.

**Alternative Tools:**

- **Bowtie2**: Another aligner that uses its own indexing method.

---

**Introduction**

Hello everyone,

Welcome back to our Next-Generation Sequencing (NGS) data analysis series. Today, we'll focus on **Step 4: Aligning Reads to the Reference Genome Using BWA MEM**. This step is pivotal as it maps your raw sequencing reads to a known reference genome, laying the foundation for all subsequent analyses, including variant calling and genomic interpretations.

We'll delve into the **BWA-MEM algorithm**, explore the **specific parameters used in the command**, and understand **why these settings are chosen** to optimize alignment accuracy and efficiency.

Let's dive in.

---

## **Step 4: Align Reads to the Reference Genome Using BWA MEM**

### **Overview**

**Read Alignment** is the process of mapping sequencing reads to a reference genome to determine their original genomic locations. Accurate alignment is crucial for reliable downstream analyses such as variant calling, structural variation detection, and gene expression profiling.

**BWA-MEM** (Burrows-Wheeler Aligner - Maximal Exact Matches) is a widely used algorithm for aligning short to moderately long sequencing reads (typically 70bp to 1Mbp) to a reference genome. It's known for its balance between speed and accuracy, making it suitable for a variety of sequencing applications.

### **BWA-MEM Algorithm**

BWA-MEM operates using the following key principles:

1. **Maximal Exact Matches (MEMs):**
   - **Definition:** Long exact matches between the read and the reference genome.
   - **Function:** Identifies seed regions in the read that exactly match regions in the reference. These seeds serve as anchors for the alignment process.
   - **Advantage:** Efficiently handles reads with high similarity to the reference, even in the presence of sequencing errors or small variations.

2. **Seed-and-Extend Strategy:**
   - **Seeding:** Finds multiple MEMs within each read.
   - **Extending:** Extends these seeds into full alignments, allowing for mismatches and gaps.
   - **Local Alignment:** Focuses on regions with high similarity, enabling accurate alignment in repetitive or complex genomic regions.

3. **Handling of Variations:**
   - **Indels and Mismatches:** Incorporates insertions, deletions, and mismatches into the alignment, accommodating structural variations and sequencing errors.
   - **Pair-End Information:** Utilizes the information from paired-end reads to improve alignment accuracy, ensuring that read pairs are correctly oriented and at expected distances.

4. **Scoring System:**
   - **Alignment Scores:** Assigns scores based on matches, mismatches, and gaps to determine the optimal alignment path.
   - **Best Alignment Selection:** Chooses the alignment with the highest score, representing the most probable origin of the read.

### **Command Breakdown and Parameters**

Let's dissect the command used for alignment:

```bash
bwa mem -M -t 6 -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tSM:sample_1' \
~/NGS_course/data/references/$ref \
~/NGS_course/data/raw/SRR26688467_1.fastq \
~/NGS_course/data/raw/SRR26688467_2.fastq \
> aligned_reads.sam
```

#### **1. `bwa mem`**

- **Function:** Invokes the BWA-MEM algorithm for aligning sequencing reads to the reference genome.
- **Usage Context:** Optimal for aligning reads ranging from 70bp to several hundred kilobases in length, commonly used for paired-end sequencing data.

#### **2. `-M`**

- **Purpose:** Marks shorter split hits as secondary alignments.
- **Detailed Explanation:**
  - **Secondary Alignments:** When a read maps to multiple locations in the genome, the best alignment is primary, and others are secondary.
  - **Compatibility:** The `-M` flag ensures compatibility with downstream tools like Picard and GATK, which expect secondary alignments to be marked appropriately.
- **Impact:** Prevents issues during downstream processing by correctly flagging multiple alignments, aiding in accurate variant calling and duplicate marking.

#### **3. `-t 6`**

- **Purpose:** Specifies the number of threads (CPU cores) to use during alignment.
- **Detailed Explanation:**
  - **Multithreading:** Enables parallel processing, significantly speeding up the alignment process, especially with large datasets.
  - **Value Selection:** `6` threads are allocated, balancing resource usage and performance. The optimal number depends on the available hardware.
- **Impact:** Reduces alignment time, enhancing computational efficiency without overburdening system resources.

#### **4. `-R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tSM:sample_1'`**

- **Purpose:** Adds read group information to each read in the SAM output.
- **Detailed Explanation:**
  - **Read Groups:** Metadata that provides context about the sequencing run, library, platform, and sample.
  - **Components:**
    - **`ID:sample_1`**: Unique identifier for the read group.
    - **`LB:sample_1`**: Library identifier, indicating the DNA library from which the reads originated.
    - **`PL:ILLUMINA`**: Platform used for sequencing, specifying the sequencing technology (e.g., Illumina).
    - **`SM:sample_1`**: Sample name, linking reads to the biological sample.
  - **Formatting:** Uses tab-separated (`\t`) fields to delineate different metadata components.
- **Impact:**
  - **Downstream Compatibility:** Essential for tools like GATK, which utilize read group information for accurate variant calling, especially in multi-sample or multiplexed sequencing scenarios.
  - **Data Organization:** Facilitates tracking of reads back to their source libraries and samples, enhancing data management and reproducibility.

#### **5. Reference Genome (`~/NGS_course/data/references/$ref`)**

- **Purpose:** Specifies the reference genome to which the reads will be aligned.
- **Detailed Explanation:**
  - **Reference Selection:** The choice of reference genome (e.g., human hg38) affects alignment accuracy and variant detection.
  - **Preparation:** Must be indexed with BWA using `bwa index` before alignment.
- **Impact:** Ensures that reads are accurately mapped to the correct genomic locations, influencing the reliability of downstream analyses.

#### **6. Input Read Files (`~/NGS_course/data/raw/SRR26688467_1.fastq` and `SRR26688467_2.fastq`)**

- **Purpose:** Provide the paired-end sequencing reads for alignment.
- **Detailed Explanation:**
  - **Paired-End Reads:** Consist of two FASTQ files representing the forward (`_1.fastq`) and reverse (`_2.fastq`) reads of each DNA fragment.
  - **Paired Information:** Enhances alignment accuracy by providing information about the insert size and orientation of read pairs.
- **Impact:** Improves the reliability of read mapping, especially in repetitive regions or areas with structural variations.

#### **7. Output Redirection (`> aligned_reads.sam`)**

- **Purpose:** Directs the alignment output to a SAM (Sequence Alignment/Map) file.
- **Detailed Explanation:**
  - **SAM Format:** A text-based format that records alignment information for each read.
  - **Usage:** Serves as an intermediate file for further processing steps like conversion to BAM, sorting, and duplicate marking.
- **Impact:** Creates a detailed record of read alignments, which is essential for quality assessment and variant calling.

### **Algorithmic Workflow of BWA-MEM in This Context**

1. **Initialization:**
   - **Reference Genome:** BWA-MEM loads the indexed reference genome (`$ref`), enabling efficient searching and alignment.
   - **Read Groups:** Integrates read group information (`-R`) into the alignment records, embedding metadata into the SAM output.

2. **Processing Paired-End Reads:**
   - **Read Pairing:** Processes forward and reverse reads as pairs, utilizing their insert size and orientation to enhance alignment accuracy.
   - **Bidirectional Alignment:** Aligns each read in both directions to identify the optimal mapping location.

3. **Identifying Maximal Exact Matches (MEMs):**
   - **Seeding:** Searches for MEMs within each read, serving as anchors for the alignment.
   - **Extension:** Extends these MEMs into full alignments, allowing for mismatches and gaps to accommodate sequencing errors and genetic variations.

4. **Scoring and Selecting Best Alignments:**
   - **Alignment Scoring:** Assigns scores based on the number of matches, mismatches, and gaps.
   - **Best Alignment Selection:** Chooses the alignment with the highest score as the primary alignment, while marking shorter or suboptimal alignments as secondary (`-M`).

5. **Output Generation:**
   - **SAM File Creation:** Writes the alignment information, including read group metadata and alignment flags, to `aligned_reads.sam`.

### **Parameter Optimization and Considerations**

1. **Thread Allocation (`-t`):**
   - **Balance:** Assigning too many threads can lead to resource contention, while too few may underutilize available CPU cores.
   - **Recommendation:** Align the number of threads with the system's CPU capacity. For example, on a machine with 8 cores, using `-t 6` provides a balance between performance and resource availability.

2. **Read Group Specification (`-R`):**
   - **Uniqueness:** Ensure that each read group has a unique `ID` to prevent conflicts in downstream analyses.
   - **Comprehensiveness:** Include all relevant metadata (e.g., library, platform, sample) to facilitate accurate variant calling and sample tracking.

3. **Handling Multiple Alignments (`-M`):**
   - **Secondary Alignments:** Properly marking secondary alignments aids tools like Picard and GATK in distinguishing between primary and alternative mappings.
   - **Impact on Variant Calling:** Prevents inflated duplicate counts and ensures accurate variant frequency estimation.

4. **Reference Genome Quality:**
   - **Completeness:** Use a high-quality, well-annotated reference genome to enhance alignment accuracy.
   - **Consistency:** Ensure that the same reference genome is used consistently across all pipeline steps to maintain data integrity.

### **Best Practices**

1. **Pre-Alignment Preparation:**
   - **Reference Indexing:** Before running BWA-MEM, index the reference genome using `bwa index` to enable efficient alignment.
     ```bash
     bwa index ~/NGS_course/data/references/hg38.fa
     ```

2. **Quality Control:**
   - **Read Quality Assessment:** Perform quality checks on raw reads (e.g., using FastQC) to identify and address issues like low-quality bases or adapter contamination before alignment.
   
3. **Resource Management:**
   - **Hardware Utilization:** Optimize thread usage based on available CPU cores to maximize alignment speed without overloading the system.
   - **Memory Allocation:** Ensure sufficient RAM is available, especially when aligning large datasets, to prevent bottlenecks.

4. **Logging and Monitoring:**
   - **Process Logging:** Use descriptive `echo` statements to log progress and facilitate troubleshooting.
   - **Error Monitoring:** Capture and review alignment logs to identify and resolve issues promptly.

5. **Reproducibility:**
   - **Script Documentation:** Clearly document all parameters and settings used in alignment commands to ensure reproducibility.
   - **Version Control:** Track the versions of tools (e.g., BWA) used in the pipeline to maintain consistency across analyses.

6. **Post-Alignment Verification:**
   - **Alignment Statistics:** Review alignment statistics (e.g., percentage of mapped reads, alignment scores) to assess alignment quality.
   - **Visualization:** Use genome browsers (e.g., IGV) to visually inspect alignments for potential anomalies or biases.

### **Example Alignment Workflow**

1. **Index the Reference Genome (if not already indexed):**
   ```bash
   bwa index ~/NGS_course/data/references/hg38.fa
   ```

2. **Run BWA-MEM Alignment:**
   ```bash
   echo "$(date "$date_format") - Aligning reads with BWA MEM"
   cd ~/NGS_course/results/alignment/
   bwa mem -M -t 6 -R '@RG\tID:sample_1\tLB:sample_1\tPL:ILLUMINA\tSM:sample_1' \
   ~/NGS_course/data/references/hg38.fa \
   ~/NGS_course/data/raw/SRR26688467_1.fastq \
   ~/NGS_course/data/raw/SRR26688467_2.fastq \
   > aligned_reads.sam
   ```

3. **Review Alignment Statistics (Optional but Recommended):**
   - **Command:**
     ```bash
     samtools flagstat aligned_reads.sam
     ```
   - **Purpose:** Provides a summary of alignment statistics, such as total reads, mapped reads, and duplicate reads.

### **Troubleshooting Tips**

1. **Low Mapping Rate:**
   - **Potential Causes:** Poor read quality, incorrect reference genome, presence of adapters or contaminants.
   - **Solutions:**
     - Perform read trimming to remove low-quality bases and adapter sequences.
     - Verify that the correct reference genome is being used.
     - Assess read quality using quality control tools like FastQC.

2. **High Number of Secondary Alignments:**
   - **Potential Causes:** Highly repetitive regions in the genome, low-quality reads.
   - **Solutions:**
     - Consider filtering out low-quality reads before alignment.
     - Use additional alignment parameters to handle repetitive sequences better.

3. **Resource Constraints:**
   - **Symptoms:** Alignment process slows down, or the system becomes unresponsive.
   - **Solutions:**
     - Reduce the number of threads (`-t`) to alleviate CPU strain.
     - Allocate more memory if possible, or process data in smaller batches.

4. **Incorrect Read Group Information:**
   - **Impact:** Downstream tools may misinterpret read origins, leading to inaccurate variant calling.
   - **Solutions:**
     - Double-check the syntax and values in the `-R` parameter.
     - Ensure that each sample or library has a unique read group identifier.

### **Conclusion**

Aligning sequencing reads to a reference genome is a foundational step in the NGS data analysis pipeline. Utilizing the **BWA-MEM algorithm** with carefully chosen parameters ensures high-quality alignments, which are crucial for accurate variant detection and subsequent analyses. By understanding the algorithmic principles and optimizing parameters like thread allocation and read group specifications, you can enhance both the efficiency and reliability of your alignment process.

Remember to adhere to best practices, such as thorough quality control, consistent reference genome usage, and diligent logging, to maintain data integrity and reproducibility throughout your analyses.

Feel free to ask any questions or request further clarification on any aspects of the BWA-MEM alignment process!

---

Here's an overview of the steps we'll cover:

1. **Step 5:** Convert SAM to Sorted BAM using Picard
2. **Step 6:** Collect Alignment and Insert Size Metrics (Optional)
3. **Step 7:** Mark Duplicates with Picard
4. **Step 8:** Create Sequence Dictionary and Index the Reference Genome
5. **Step 9:** Index the Known Sites VCF File

Let's begin.

---

## **Step 5: Convert SAM to Sorted BAM Using Picard**

### **Overview**

After aligning your sequencing reads to the reference genome using tools like BWA-MEM (as in Step 4), you obtain a SAM (Sequence Alignment/Map) file. SAM files are text-based and can be quite large, making them inefficient for storage and processing. To address this, we convert the SAM file to a BAM (Binary Alignment/Map) file, which is the compressed binary version of SAM. Additionally, we sort the BAM file by genomic coordinates to facilitate efficient data retrieval and compatibility with downstream tools.

### **Command**

```bash
echo "$(date "$date_format") - Sorting SAM file and converting to BAM"
$picard SortSam \
    INPUT=aligned_reads.sam \
    OUTPUT=sorted_reads.bam \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT
```

### **Breaking Down the Command**

- **`echo "$(date "$date_format") - Sorting SAM file and converting to BAM"`**
  - **Purpose:** Logs the current date and a descriptive message indicating the start of the SAM to BAM conversion process.
  
- **`$picard SortSam`**
  - **Description:** Invokes the `SortSam` tool from the Picard suite, which sorts and converts SAM files to BAM format.
  - **`$picard`** is typically an environment variable pointing to the Picard executable or a wrapper script.

- **Parameters:**
  - **`INPUT=aligned_reads.sam`**
    - **Description:** Specifies the input SAM file generated from the alignment step.
  
  - **`OUTPUT=sorted_reads.bam`**
    - **Description:** Defines the name of the output BAM file that will contain the sorted alignments.
  
  - **`SORT_ORDER=coordinate`**
    - **Description:** Instructs Picard to sort the alignments based on their genomic coordinates (i.e., by chromosome and position).
  
  - **`VALIDATION_STRINGENCY=SILENT`**
    - **Description:** Sets the validation stringency to `SILENT`, meaning Picard will not report warnings about minor formatting issues in the SAM file.
    - **Alternative Values:** `STRICT`, `LENIENT`, etc., which control how strictly the tool validates the input files.

### **Associated File Formats**

1. **SAM File (`.sam`):**
   - **Description:** Text-based format for storing aligned sequencing reads.
   - **Structure:** Each alignment is represented by a single line with multiple tab-separated fields.
   - **Example Entry:**
     ```
     SRR26688467.1    99    chr1    10000    60    76M    =    10100    200    ATGCGTACGTTAGCTAGCTAGCTAGCTAGCTA    BBBFFFFFFFFFFFFFFFFFFFFFFFFFFFF    AS:i:0    XN:i:0
     ```

2. **BAM File (`.bam`):**
   - **Description:** Binary, compressed version of the SAM format.
   - **Advantages:**
     - **Storage Efficiency:** Significantly smaller in size compared to SAM.
     - **Speed:** Faster processing and data retrieval.
     - **Compatibility:** Required by many downstream tools (e.g., Picard, GATK).
   - **Example Usage:**
     - Viewing a BAM file in a human-readable format:
       ```bash
       samtools view sorted_reads.bam | head
       ```

### **Why This Step Is Important**

1. **Storage Efficiency:**
   - SAM files are large and cumbersome to work with. Converting to BAM reduces file size, saving storage space and enabling faster data processing.

2. **Sorting Alignments:**
   - Sorting alignments by genomic coordinates is essential for many downstream analyses, including duplicate marking, variant calling, and visualization in genome browsers.

3. **Compatibility:**
   - Many bioinformatics tools require sorted BAM files as input. Ensuring your data is in the correct format streamlines the pipeline and avoids compatibility issues.

### **Best Practices**

- **Consistent Naming Conventions:**
  - Use clear and consistent file naming to track data through the pipeline (e.g., `sorted_reads.bam`).

- **Resource Allocation:**
  - Ensure sufficient memory and CPU resources are available, especially when working with large SAM files.

- **Validation Stringency:**
  - While `SILENT` reduces log verbosity, consider using `LENIENT` or `STRICT` during initial runs to catch potential issues in the SAM file.

- **Indexing:**
  - After sorting, consider indexing the BAM file using `samtools index` for faster access in downstream analyses.

### **Example Workflow**

1. **Sort and Convert SAM to BAM:**
   ```bash
   picard SortSam \
       INPUT=aligned_reads.sam \
       OUTPUT=sorted_reads.bam \
       SORT_ORDER=coordinate \
       VALIDATION_STRINGENCY=SILENT
   ```

2. **Index the Sorted BAM File:**
   ```bash
   samtools index sorted_reads.bam
   ```

---

## **Step 6: Collect Alignment and Insert Size Metrics (Optional)**

### **Overview**

Collecting alignment and insert size metrics provides valuable insights into the quality and characteristics of your sequencing data. These metrics help assess the efficiency of the alignment process, the quality of the sequencing library, and potential biases that may affect downstream analyses like variant calling.

### **Commands**

1. **Collect Alignment Summary Metrics:**

   ```bash
   echo "$(date "$date_format") - Collecting alignment metrics"
   $picard CollectAlignmentSummaryMetrics \
       R=~/NGS_course/data/references/$ref \
       I=sorted_reads.bam \
       O=alignment_metrics.txt
   ```

2. **Collect Insert Size Metrics:**

   ```bash
   echo "$(date "$date_format") - Collecting insert size metrics"
   $picard CollectInsertSizeMetrics \
       INPUT=sorted_reads.bam \
       OUTPUT=insert_metrics.txt \
       HISTOGRAM_FILE=insert_size_histogram.pdf
   ```

### **Breaking Down the Commands**

#### **1. CollectAlignmentSummaryMetrics**

- **`$picard CollectAlignmentSummaryMetrics`**
  - **Description:** Picard tool that gathers various alignment statistics from a BAM file.
  
- **Parameters:**
  - **`R=~/NGS_course/data/references/$ref`**
    - **Description:** Reference genome file used for alignment.
  
  - **`I=sorted_reads.bam`**
    - **Description:** Input sorted BAM file containing aligned reads.
  
  - **`O=alignment_metrics.txt`**
    - **Description:** Output text file that will store the collected alignment metrics.

#### **2. CollectInsertSizeMetrics**

- **`$picard CollectInsertSizeMetrics`**
  - **Description:** Picard tool that computes insert size metrics from a BAM file, providing insights into the sequencing library's fragment size distribution.
  
- **Parameters:**
  - **`INPUT=sorted_reads.bam`**
    - **Description:** Input sorted BAM file.
  
  - **`OUTPUT=insert_metrics.txt`**
    - **Description:** Output text file containing numerical insert size metrics.
  
  - **`HISTOGRAM_FILE=insert_size_histogram.pdf`**
    - **Description:** Output PDF file that visualizes the insert size distribution as a histogram.

### **Associated File Formats**

1. **Metrics Files (`.txt`):**
   - **Description:** Text files containing tabulated statistics.
   - **Examples:**
     - **Alignment Metrics (`alignment_metrics.txt`):**
       ```
       CATEGORY	TOTAL_READS	PF_READS	UNPAIRED_READS	...
       UNPAIRED_READS	1000000	990000	...
       PAIR_READS	2000000	1980000	...
       ```
     - **Insert Size Metrics (`insert_metrics.txt`):**
       ```
       MEAN_INSERT_SIZE	INSERT_SIZE_STDDEV	...
       300	50	...
       ```

2. **Histogram Files (`.pdf`):**
   - **Description:** Graphical representations of data distributions.
   - **Example:**
     - **Insert Size Histogram (`insert_size_histogram.pdf`):**
       - A bar graph showing the frequency of different insert sizes, allowing visual assessment of the library's fragment size distribution.

### **Why This Step Is Important**

1. **Quality Assessment:**
   - **Alignment Metrics:** Provide information on the proportion of reads successfully aligned, mapping quality distributions, and potential biases.
   - **Insert Size Metrics:** Reveal the consistency of fragment sizes, which is crucial for paired-end sequencing and accurate variant calling.

2. **Library Evaluation:**
   - **Fragment Size Distribution:** Ensures that the library preparation was successful and that the insert sizes are within expected ranges.

3. **Troubleshooting:**
   - **Identify Issues:** Metrics can help identify problems like excessive duplicate reads, uneven coverage, or unexpected insert size distributions that may require adjustments in the protocol.

### **Understanding the Metrics**

#### **1. Alignment Summary Metrics**

Key metrics typically include:

- **Total Reads:** The total number of reads processed.
- **PF Reads:** Reads that pass the vendor filter.
- **Mapped Reads:** Number and percentage of reads that successfully aligned to the reference genome.
- **Properly Paired:** Percentage of read pairs where both reads are properly aligned.
- **Mismatches per Read:** Average number of mismatches per read.
- **Error Rate:** Overall error rate in the alignment.

#### **2. Insert Size Metrics**

Key metrics typically include:

- **Mean Insert Size:** The average distance between the start of the first read and the end of the second read in a pair.
- **Insert Size Standard Deviation:** Measures the variability in insert sizes.
- **Median Insert Size:** The middle value of the insert size distribution.
- **Histogram Visualization:** Provides a visual representation of insert size distribution, helping identify bimodal distributions or unexpected peaks.

### **Best Practices**

1. **Regular Monitoring:**
   - Collect metrics at various stages of the pipeline to monitor data quality continuously.

2. **Consistent Reference Files:**
   - Ensure that the reference genome used for collecting metrics matches the one used for alignment to avoid discrepancies.

3. **Interpretation of Metrics:**
   - Compare collected metrics against expected standards or benchmarks to assess data quality.
   - For insert size metrics, expect a unimodal distribution centered around the library's target fragment size.

4. **Documentation:**
   - Keep records of metrics for reproducibility and future reference, especially when troubleshooting or optimizing protocols.

### **Example Workflow**

1. **Collect Alignment Summary Metrics:**
   ```bash
   picard CollectAlignmentSummaryMetrics \
       R=~/NGS_course/data/references/hg38.fa \
       I=sorted_reads.bam \
       O=alignment_metrics.txt
   ```

2. **Collect Insert Size Metrics:**
   ```bash
   picard CollectInsertSizeMetrics \
       INPUT=sorted_reads.bam \
       OUTPUT=insert_metrics.txt \
       HISTOGRAM_FILE=insert_size_histogram.pdf
   ```

3. **Review Metrics:**
   - Open `alignment_metrics.txt` and `insert_metrics.txt` to assess alignment efficiency and insert size distribution.
   - View `insert_size_histogram.pdf` to visualize the insert size distribution.

---

## **Step 7: Mark Duplicates with Picard**

### **Overview**

During the sequencing process, especially during PCR amplification, duplicate reads can be introduced. These duplicates are multiple copies of the same original DNA fragment and do not represent independent observations. If not addressed, duplicates can bias variant calling by artificially inflating allele counts, leading to false-positive variant calls. **Marking duplicates** is essential to mitigate this bias without removing the reads, allowing downstream tools to account for duplicates appropriately.

### **Command**

```bash
echo "$(date "$date_format") - Marking duplicates"
$picard MarkDuplicates \
    INPUT=sorted_reads.bam \
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT
```

### **Breaking Down the Command**

- **`echo "$(date "$date_format") - Marking duplicates"`**
  - **Purpose:** Logs the current date and a descriptive message indicating the start of the duplicate marking process.

- **`$picard MarkDuplicates`**
  - **Description:** Invokes the `MarkDuplicates` tool from the Picard suite, which identifies and marks duplicate reads in a BAM file.

- **Parameters:**
  - **`INPUT=sorted_reads.bam`**
    - **Description:** Specifies the input BAM file that has been sorted by genomic coordinates.
  
  - **`OUTPUT=dedup_reads.bam`**
    - **Description:** Defines the name of the output BAM file with duplicates marked.
  
  - **`METRICS_FILE=metrics.txt`**
    - **Description:** Output file that will contain statistics about the duplication process, including the number of duplicates found.
  
  - **`CREATE_INDEX=true`**
    - **Description:** Instructs Picard to create an index file (`.bai`) for the output BAM file.
  
  - **`VALIDATION_STRINGENCY=SILENT`**
    - **Description:** Sets the validation stringency to `SILENT`, suppressing warnings about minor formatting issues in the BAM file.

### **Associated File Formats**

1. **Deduplicated BAM File (`dedup_reads.bam`):**
   - **Description:** A BAM file identical to the input but with duplicate reads marked.
   - **Marking Method:** Duplicates are flagged using the `0x400` bitwise flag in the BAM format, indicating that the read is a duplicate.
   - **Usage:** Input for downstream tools like GATK HaplotypeCaller, which can ignore or appropriately handle marked duplicates during variant calling.

2. **Metrics File (`metrics.txt`):**
   - **Description:** A text file containing detailed statistics about the duplication process.
   - **Content:**
     - **Total Reads:** Total number of reads processed.
     - **Duplicate Reads:** Number and percentage of reads identified as duplicates.
     - **Unique Reads:** Number of reads considered unique.
   - **Example Content:**
     ```
     METRIC\tTYPE\tPROGRAM\tRUNTIME_SECONDS\tCATEGORY\tSTART\tSTOP\tVALUE
     METRIC\tREADS\tMarkDuplicates\t120\tTOTAL_READS\t0\t120\t2000000
     METRIC\tREADS\tMarkDuplicates\t120\tUNIQUE_READS\t0\t120\t1800000
     METRIC\tREADS\tMarkDuplicates\t120\tDUPLICATE_READS\t0\t120\t200000
     ```

3. **BAM Index File (`dedup_reads.bam.bai`):**
   - **Description:** Binary index file that allows rapid access to specific regions within the BAM file.
   - **Usage:** Required by tools that process BAM files to enable efficient data retrieval.

### **Why This Step Is Important**

1. **Preventing Bias in Variant Calling:**
   - Duplicate reads can falsely increase the confidence in allele counts, leading to the incorrect identification of variants.
   - Marking duplicates allows variant callers to recognize and handle these reads appropriately, reducing false-positive variant calls.

2. **Maintaining Data Integrity:**
   - Instead of removing duplicates, marking preserves all reads while providing information about their duplication status.
   - This approach ensures that no potential data is lost, allowing for flexibility in downstream analyses.

3. **Enhancing Downstream Analyses:**
   - Many downstream tools are designed to recognize and ignore marked duplicates, ensuring accurate results.
   - Proper duplicate marking is a prerequisite for high-quality variant calling and other analyses.

### **Understanding Duplicate Marking**

- **Identification of Duplicates:**
  - Reads are considered duplicates if they originate from the same original DNA fragment, typically having identical start and end positions in the genome.
  
- **Marking vs. Removal:**
  - **Marking:** Flags duplicate reads without removing them, preserving all data.
  - **Removal:** Physically deletes duplicate reads, which is generally not recommended as it may remove useful information.
  
- **Flags in BAM:**
  - **`0x400` Flag:** Indicates that a read is a duplicate.

### **Best Practices**

1. **Consistent Read Group Information:**
   - Ensure that read groups are correctly defined during alignment, as `MarkDuplicates` relies on read group information to accurately identify duplicates across different libraries or sequencing runs.

2. **Resource Allocation:**
   - Duplicate marking can be resource-intensive for large BAM files. Allocate sufficient memory and CPU resources to handle the process efficiently.

3. **Review Metrics:**
   - After running `MarkDuplicates`, review the `metrics.txt` file to assess the duplication rate. High duplication rates may indicate issues with library preparation, such as excessive PCR amplification.

4. **Avoid Removing Duplicates:**
   - Instead of removing duplicates, prefer marking them. This approach maintains data integrity and allows downstream tools to handle duplicates appropriately.

### **Example Workflow**

1. **Mark Duplicates:**
   ```bash
   picard MarkDuplicates \
       INPUT=sorted_reads.bam \
       OUTPUT=dedup_reads.bam \
       METRICS_FILE=metrics.txt \
       CREATE_INDEX=true \
       VALIDATION_STRINGENCY=SILENT
   ```

2. **Review Duplication Metrics:**
   - Open and examine `metrics.txt` to understand the extent of duplication.
     ```bash
     cat metrics.txt
     ```
     - **Interpretation:** Look for the percentage of duplicate reads. A high percentage (e.g., >20%) may warrant further investigation into library preparation protocols.

3. **Index the Deduplicated BAM File:**
   - If `CREATE_INDEX=true` was used, the index file (`dedup_reads.bam.bai`) is created automatically.
   - If not, index manually:
     ```bash
     samtools index dedup_reads.bam
     ```

---

## **Step 8: Create Sequence Dictionary and Index the Reference Genome**

### **Overview**

Creating a **sequence dictionary** and indexing the reference genome are essential steps to prepare the reference for efficient access during variant calling and other analyses. These preparatory steps ensure that bioinformatics tools can quickly and accurately reference the genome during data processing.

### **Commands**

1. **Create Sequence Dictionary:**

   ```bash
   echo "$(date "$date_format") - Creating sequence dictionary"
   $picard CreateSequenceDictionary \
       R=~/NGS_course/data/references/$ref \
       O=~/NGS_course/data/references/${ref%.fa}.dict
   ```

2. **Index the Reference Genome with Samtools:**

   ```bash
   echo "$(date "$date_format") - Indexing the reference genome with samtools"
   samtools faidx ~/NGS_course/data/references/$ref
   ```

### **Breaking Down the Commands**

#### **1. CreateSequenceDictionary**

- **`$picard CreateSequenceDictionary`**
  - **Description:** Picard tool that generates a sequence dictionary for the reference genome.
  
- **Parameters:**
  - **`R=~/NGS_course/data/references/$ref`**
    - **Description:** Specifies the input reference genome FASTA file.
  
  - **`O=~/NGS_course/data/references/${ref%.fa}.dict`**
    - **Description:** Defines the output sequence dictionary file, typically with a `.dict` extension.
    - **`${ref%.fa}`**: This is a shell parameter expansion that removes the `.fa` extension from the reference file name.

#### **2. samtools faidx**

- **`samtools faidx ~/NGS_course/data/references/$ref`**
  - **Description:** Uses Samtools to index the reference genome FASTA file.
  
- **Parameters:**
  - **`~/NGS_course/data/references/$ref`**
    - **Description:** Specifies the reference genome FASTA file to be indexed.
  
- **Output:**
  - **`.fai` File:** Samtools creates an index file with the `.fai` extension, which maps each reference sequence name to its position in the FASTA file.

### **Associated File Formats**

1. **Sequence Dictionary (`.dict`):**
   - **Description:** A text file that provides a mapping of sequence names and lengths in the reference genome.
   - **Structure:**
     - **Header Line:**
       ```
       @HD VN:1.6 SO:coordinate
       ```
     - **Sequence Lines:**
       ```
       @SQ SN:chr1 LN:248956422 UR:file:/path/to/reference.fa
       @SQ SN:chr2 LN:242193529 UR:file:/path/to/reference.fa
       ```
       - **`@SQ`:** Indicates a sequence dictionary entry.
       - **`SN:`**: Sequence name (e.g., chromosome name).
       - **`LN:`**: Length of the sequence.
       - **`UR:`**: URI to the reference FASTA file (optional but recommended).
   
   - **Example Entry:**
     ```
     @SQ	SN:chr1	LN:248956422	UR:file:/data/references/hg38.fa
     ```

2. **FASTA Index File (`.fai`):**
   - **Description:** An index file created by Samtools that enables rapid access to specific regions of the reference genome without reading the entire file.
   - **Structure:** Tab-separated values, one line per reference sequence.
     ```
     chr1	248956422	52	60	61
     chr2	242193529	253404903	60	61
     ```
     - **Fields:**
       1. **Sequence Name:** `chr1`
       2. **Sequence Length:** `248956422`
       3. **Offset:** Byte offset in the FASTA file where the sequence starts.
       4. **Line Bases:** Number of bases per line in the FASTA file.
       5. **Line Width:** Number of bytes per line (including newline characters).
   
   - **Example Entry:**
     ```
     chr1	248956422	52	60	61
     ```

### **Why This Step Is Important**

1. **Efficient Data Access:**
   - **Sequence Dictionary:** Facilitates tools like GATK to understand the structure and naming of sequences in the reference genome.
   - **FASTA Index:** Allows tools to quickly retrieve specific genomic regions without loading the entire reference genome into memory, significantly speeding up analyses.

2. **Tool Compatibility:**
   - Many bioinformatics tools, including GATK and Picard, require a sequence dictionary and indexed reference genome to function correctly.

3. **Data Integrity:**
   - Ensures that all tools reference the same version and build of the genome, preventing inconsistencies and errors in analyses.

### **Understanding the File Formats**

#### **1. Sequence Dictionary (`.dict`):**

- **Purpose:** Provides metadata about each sequence (e.g., chromosomes) in the reference genome.
- **Usage:** Required by GATK for variant calling and other analyses to map read alignments correctly.
- **Key Information:**
  - **Sequence Names:** Must match exactly between the sequence dictionary and the reference FASTA file.
  - **Sequence Lengths:** Accurate lengths ensure correct positioning of reads and variants.

#### **2. FASTA Index (`.fai`):**

- **Purpose:** Enables random access to any region of the reference genome.
- **Usage:** Used by Samtools and other tools to fetch specific genomic regions efficiently.
- **Key Information:**
  - **Sequence Name and Length:** Aligns with the reference FASTA file.
  - **Offset and Line Information:** Critical for tools to calculate byte positions and retrieve data accurately.

### **Best Practices**

1. **Consistency:**
   - Ensure that the reference genome FASTA file used for creating the sequence dictionary and index matches exactly with the one used in the alignment step.
  
2. **Version Control:**
   - Keep track of reference genome versions to maintain reproducibility. Different builds (e.g., hg38 vs. hg19) can lead to discrepancies in analyses.
  
3. **File Management:**
   - Store the sequence dictionary and index files in the same directory as the reference FASTA file for easy access and organization.
  
4. **Validation:**
   - Use tools like GATK's `ValidateSamFile` to ensure that the sequence dictionary is correctly formatted and matches the reference genome.

### **Example Workflow**

1. **Create Sequence Dictionary:**
   ```bash
   picard CreateSequenceDictionary \
       R=~/NGS_course/data/references/hg38.fa \
       O=~/NGS_course/data/references/hg38.dict
   ```

2. **Index Reference Genome with Samtools:**
   ```bash
   samtools faidx ~/NGS_course/data/references/hg38.fa
   ```

3. **Verify Files:**
   - **Sequence Dictionary (`hg38.dict`):**
     ```
     @HD	VN:1.6	SO:coordinate
     @SQ	SN:chr1	LN:248956422	UR:file:/data/references/hg38.fa
     @SQ	SN:chr2	LN:242193529	UR:file:/data/references/hg38.fa
     ...
     ```
   
   - **FASTA Index (`hg38.fa.fai`):**
     ```
     chr1	248956422	52	60	61
     chr2	242193529	253404903	60	61
     ...
     ```

---

## **Step 9: Index the Known Sites VCF File**

### **Overview**

Indexing the **Known Sites VCF File** is a preparatory step for processes like Base Quality Score Recalibration (BQSR) and variant calling. Known sites (e.g., from databases like dbSNP or Mills) provide a list of known genetic variants that help distinguish true variants from sequencing errors. Indexing the VCF file ensures efficient access to specific genomic regions during analysis, which is crucial for performance and accuracy.

### **Command**

```bash
echo "$(date "$date_format") - Indexing known sites VCF file"
cd ~/NGS_course/data/references/
if [ ! -f "${knownsites}.tbi" ]; then
    tabix -p vcf $knownsites
fi
```

### **Breaking Down the Command**

- **`echo "$(date "$date_format") - Indexing known sites VCF file"`**
  - **Purpose:** Logs the current date and a descriptive message indicating the start of the VCF indexing process.

- **`cd ~/NGS_course/data/references/`**
  - **Purpose:** Changes the working directory to where the known sites VCF file is located.

- **`if [ ! -f "${knownsites}.tbi" ]; then`**
  - **Description:** Checks if the index file (`.tbi`) for the known sites VCF already exists.
  
- **`tabix -p vcf $knownsites`**
  - **Description:** Uses Tabix to index the VCF file.
  
- **`fi`**
  - **Description:** Ends the conditional statement.

### **Associated File Formats**

1. **VCF File (`.vcf` or `.vcf.gz`):**
   - **Description:** Variant Call Format file containing known genetic variants.
   - **Usage:** Input for BQSR and variant calling tools to provide a reference of known variants.
   - **Example Entry:**
     ```
     #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
     chr1	10177	rs367896724	A	AC	100	PASS	AC=1;AF=0.5;AN=2
     ```

2. **VCF Index File (`.tbi`):**
   - **Description:** Tabix index file for compressed VCF files (`.vcf.gz`).
   - **Usage:** Enables rapid querying of specific genomic regions within the VCF file.
   - **Example Entry:**
     - The `.tbi` file is a binary file and does not contain human-readable content. It is used by tools like Tabix and GATK to access variant information efficiently.

### **Why This Step Is Important**

1. **Efficiency in Data Access:**
   - **Indexed VCF Files:** Allow tools to quickly retrieve variants within specific genomic regions without scanning the entire file, significantly speeding up analyses.
   
2. **Essential for BQSR and Variant Calling:**
   - **BQSR:** Uses known variants to identify and correct systematic errors in base quality scores.
   - **Variant Calling:** References known variants to differentiate between true genetic variants and sequencing artifacts.

3. **Resource Optimization:**
   - **Performance:** Indexing minimizes computational load and memory usage by enabling targeted data access.

### **Understanding Tabix and VCF Indexing**

- **Tabix:**
  - **Description:** A tool for indexing and querying BGZF-compressed, tab-delimited files like VCF.
  - **Functionality:** Creates an index file (`.tbi`) that maps genomic regions to byte offsets in the compressed VCF file.
  
- **VCF Indexing:**
  - **Requirements:**
    - **Compressed VCF:** The VCF file must be compressed using BGZF (Block GZIP), typically resulting in a `.vcf.gz` file.
    - **Sorted VCF:** The VCF file must be sorted by chromosome and position to create an effective index.
  
  - **Process:**
    1. **Compress the VCF File:**
       ```bash
       bgzip known_sites.vcf
       ```
       - **Result:** `known_sites.vcf.gz`
    
    2. **Sort the VCF File:**
       ```bash
       bcftools sort known_sites.vcf.gz -o known_sites.sorted.vcf.gz
       ```
       - **Result:** `known_sites.sorted.vcf.gz`
    
    3. **Index the Sorted VCF File:**
       ```bash
       tabix -p vcf known_sites.sorted.vcf.gz
       ```
       - **Result:** `known_sites.sorted.vcf.gz.tbi`

### **Best Practices**

1. **Ensure VCF Compression:**
   - Use BGZF compression (`.vcf.gz`) before indexing with Tabix. Regular GZIP compression (`.gz`) is not compatible with Tabix.
   - **Command:**
     ```bash
     bgzip known_sites.vcf
     ```

2. **Sort the VCF File:**
   - Before indexing, sort the VCF file by chromosome and position to comply with Tabix requirements.
   - **Command:**
     ```bash
     bcftools sort known_sites.vcf.gz -o known_sites.sorted.vcf.gz
     ```

3. **Verify the Index:**
   - After indexing, ensure that the `.tbi` file is present and correctly associated with the `.vcf.gz` file.
   - **Command:**
     ```bash
     ls known_sites.sorted.vcf.gz.tbi
     ```

4. **Automate the Process:**
   - Incorporate checks in your scripts to prevent redundant indexing, as shown in the provided command block (`if [ ! -f "${knownsites}.tbi" ]; then ... fi`).

5. **Maintain Consistency:**
   - Use the same reference genome build and known sites VCF across all samples to ensure consistency in analyses.

6. **Documentation:**
   - Keep records of the versions and sources of known sites databases used, as these can impact variant calling accuracy.

### **Example Workflow**

1. **Compress the Known Sites VCF File:**
   ```bash
   bgzip known_sites.vcf
   ```
   - **Result:** `known_sites.vcf.gz`

2. **Sort the Compressed VCF File:**
   ```bash
   bcftools sort known_sites.vcf.gz -o known_sites.sorted.vcf.gz
   ```
   - **Result:** `known_sites.sorted.vcf.gz`

3. **Index the Sorted VCF File:**
   ```bash
   tabix -p vcf known_sites.sorted.vcf.gz
   ```
   - **Result:** `known_sites.sorted.vcf.gz.tbi`

4. **Integrate into Pipeline:**
   - Use the sorted and indexed VCF file in subsequent steps like BQSR.
   - **Example Command:**
     ```bash
     gatk BaseRecalibrator \
         -I dedup_reads.bam \
         -R hg38.fa \
         --known-sites known_sites.sorted.vcf.gz \
         -O recal_data.table
     ```

### **Troubleshooting Tips**

1. **Common Errors:**
   - **Unsorted VCF File:**
     - **Error Message:** "ERROR: Could not determine the start of the index region."
     - **Solution:** Ensure the VCF file is sorted by chromosome and position before indexing.
   
   - **Incorrect Compression:**
     - **Error Message:** "ERROR: Unknown file format."
     - **Solution:** Verify that the VCF file is BGZF-compressed (`.vcf.gz`) and not regular GZIP-compressed (`.gz`).

2. **Validation:**
   - Use tools like `bcftools` to validate the integrity of the compressed and sorted VCF file.
     ```bash
     bcftools view known_sites.sorted.vcf.gz | head
     ```

3. **Permission Issues:**
   - Ensure you have the necessary read/write permissions in the directory where indexing is performed.

### **Conclusion of Step 9**

Indexing the known sites VCF file is a straightforward yet essential step that enhances the efficiency and accuracy of downstream analyses like Base Quality Score Recalibration and variant calling. By ensuring that the known sites are properly compressed, sorted, and indexed, you facilitate rapid access to variant information, thereby optimizing your NGS data analysis pipeline.

---

## **Summary of Steps 5 to 9**

| **Step** | **Description** | **Input File(s)** | **Output File(s)** |
|----------|-----------------|--------------------|---------------------|
| **5**    | Convert SAM to sorted BAM using Picard | `aligned_reads.sam` | `sorted_reads.bam` |
| **6**    | Collect alignment and insert size metrics | `sorted_reads.bam` | `alignment_metrics.txt`, `insert_metrics.txt`, `insert_size_histogram.pdf` |
| **7**    | Mark duplicates with Picard | `sorted_reads.bam` | `dedup_reads.bam`, `metrics.txt` |
| **8**    | Create sequence dictionary and index the reference genome | `hg38.fa` | `hg38.dict`, `hg38.fa.fai` |
| **9**    | Index the known sites VCF file | `known_sites.sorted.vcf.gz` | `known_sites.sorted.vcf.gz.tbi` |

---

## **Key Takeaways**

1. **Efficient Data Management:**
   - Converting SAM to BAM and sorting by coordinates enhances storage efficiency and compatibility with downstream tools.

2. **Quality Assessment:**
   - Collecting alignment and insert size metrics provides insights into data quality, aiding in the identification of potential issues early in the pipeline.

3. **Duplication Handling:**
   - Marking duplicates prevents bias in variant calling without losing valuable read data.

4. **Reference Genome Preparation:**
   - Creating a sequence dictionary and indexing the reference genome are prerequisites for accurate and efficient data processing in variant calling and other analyses.

5. **Known Sites Indexing:**
   - Indexing known variant sites ensures efficient access during Base Quality Score Recalibration and improves the accuracy of variant calls by distinguishing true variants from sequencing errors.

---

## **Best Practices Across Steps 5 to 9**

- **Consistent Reference Files:**
  - Use the same reference genome version across all steps to maintain consistency and accuracy.
  
- **Resource Management:**
  - Allocate sufficient computational resources (CPU, memory) to handle large BAM and VCF files efficiently.
  
- **Automation and Scripting:**
  - Automate repetitive tasks using shell scripts or workflow managers (e.g., Snakemake, Nextflow) to ensure reproducibility and reduce manual errors.
  
- **Documentation:**
  - Keep detailed records of all commands, parameters, and file versions used in the pipeline for transparency and reproducibility.
  
- **Regular Monitoring:**
  - Continuously monitor metrics and logs to detect and address issues promptly.

- **Backup Critical Files:**
  - Maintain backups of raw and intermediate files to prevent data loss and facilitate troubleshooting.

---

## **Questions and Discussion**

Feel free to ask any questions you might have about Steps 5 to 9, the commands and parameters used, the associated file formats, or best practices for optimizing your NGS data analysis pipeline. Understanding these steps in detail is fundamental to conducting accurate and reliable genomic analyses.

---

## **Additional Resources**

- **Picard Tools Documentation:**
  - [SortSam](https://broadinstitute.github.io/picard/command-line-overview.html#SortSam)
  - [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)
  - [CreateSequenceDictionary](https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)

- **Samtools Documentation:**
  - [Samtools faidx](http://www.htslib.org/doc/samtools-faidx.html)
  - [Samtools index](http://www.htslib.org/doc/samtools-index.html)

- **Tabix Documentation:**
  - [Tabix Usage](https://www.htslib.org/doc/tabix.html)

- **GATK Best Practices:**
  - [GATK Best Practices Pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-GATK-Best-Practices-Workflows)

- **VCF Format Specification:**
  - [VCF v4.2 Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

- **Tutorials and Workshops:**
  - [GATK Tutorials](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-GATK-Tutorials)
  - [Picard Tutorials](https://broadinstitute.github.io/picard/)



We'll delve into two pivotal steps of the NGS pipeline:

1. **Step 10: Base Quality Score Recalibration (BQSR) with GATK**
2. **Step 11: Variant Calling with GATK HaplotypeCaller**

Understanding these steps is crucial for ensuring the accuracy and reliability of your variant calls, which are fundamental for downstream analyses like association studies and functional genomics. We'll explore the commands involved, the purpose of each step, the associated file formats, and best practices to optimize your analysis.

---

## **Step 10: Base Quality Score Recalibration (BQSR) with GATK**

### **Overview**

**Base Quality Score Recalibration (BQSR)** is a crucial preprocessing step in the NGS data analysis pipeline. It aims to correct systematic errors made by the sequencing machine when it assigns quality scores to each base call. Accurate base quality scores are essential for reliable variant calling, as they influence the confidence in detecting true genetic variants versus sequencing errors.

### **Why BQSR Is Important**

1. **Corrects Systematic Errors:**
   - Sequencing machines can introduce biases in base quality scores based on various factors like machine cycle, sequence context, and read position.
   - These biases can lead to inaccurate variant calls if not corrected.

2. **Improves Variant Calling Accuracy:**
   - Accurate quality scores help variant callers distinguish between true variants and sequencing errors, reducing false positives and negatives.

3. **Enhances Downstream Analyses:**
   - Reliable variant calls are foundational for all subsequent analyses, including genotype-phenotype associations and functional studies.

### **Associated File Formats**

- **BAM (`.bam`):**
  - **Description:** Binary version of SAM files containing aligned sequencing reads.
  - **Usage:** Input for BQSR tools.
  
- **Recalibration Table (`.table`):**
  - **Description:** A text file containing recalibration data generated by `BaseRecalibrator`.
  - **Usage:** Input for `ApplyBQSR` to adjust base quality scores.
  
- **Recalibrated BAM (`.bam`):**
  - **Description:** BAM file with adjusted base quality scores.
  - **Usage:** Input for variant calling tools like HaplotypeCaller.

### **Commands and Parameters**

Let's examine the commands used for BQSR in detail.

#### **Command 1: BaseRecalibrator**

```bash
gatk BaseRecalibrator \
    -I dedup_reads.bam \
    -R ~/NGS_course/data/references/$ref \
    --known-sites $knownsites \
    -O recal_data.table
```

**Breaking Down the Command:**

- **`gatk BaseRecalibrator`:**
  - Invokes the `BaseRecalibrator` tool from the Genome Analysis Toolkit (GATK).

- **Parameters:**
  - **`-I dedup_reads.bam`**
    - **Description:** Specifies the input BAM file with duplicates marked.
    - **Importance:** Ensures that duplicate reads do not skew recalibration metrics.
  
  - **`-R ~/NGS_course/data/references/$ref`**
    - **Description:** Specifies the reference genome file.
    - **Importance:** Essential for aligning reads and identifying variants accurately.
  
  - **`--known-sites $knownsites`**
    - **Description:** Provides a VCF file of known variant sites (e.g., dbSNP, Mills).
    - **Importance:** Helps distinguish true variants from sequencing errors by providing known polymorphic sites.
  
  - **`-O recal_data.table`**
    - **Description:** Defines the output file for the recalibration table.
    - **Importance:** This table contains statistical models of the systematic errors, which are used in the next step.

**Example:**

```bash
gatk BaseRecalibrator \
    -I dedup_reads.bam \
    -R ~/NGS_course/data/references/hg38.fa \
    --known-sites ~/NGS_course/data/references/dbsnp_138.hg38.vcf.gz \
    -O recal_data.table
```

#### **Command 2: ApplyBQSR**

```bash
gatk ApplyBQSR \
    -R ~/NGS_course/data/references/$ref \
    -I dedup_reads.bam \
    --bqsr-recal-file recal_data.table \
    -O recal_reads.bam
```

**Breaking Down the Command:**

- **`gatk ApplyBQSR`:**
  - Invokes the `ApplyBQSR` tool from GATK.

- **Parameters:**
  - **`-R ~/NGS_course/data/references/$ref`**
    - **Description:** Specifies the reference genome file.
  
  - **`-I dedup_reads.bam`**
    - **Description:** Input BAM file with duplicates marked.
  
  - **`--bqsr-recal-file recal_data.table`**
    - **Description:** Provides the recalibration table generated by `BaseRecalibrator`.
  
  - **`-O recal_reads.bam`**
    - **Description:** Defines the output BAM file with recalibrated base quality scores.

**Example:**

```bash
gatk ApplyBQSR \
    -R ~/NGS_course/data/references/hg38.fa \
    -I dedup_reads.bam \
    --bqsr-recal-file recal_data.table \
    -O recal_reads.bam
```

### **Step-by-Step Explanation**

1. **BaseRecalibrator:**
   - **Function:** Analyzes patterns of covariation between the recorded quality scores and the empirical error rate.
   - **Process:**
     - Examines the BAM file to identify systematic errors in base quality scores.
     - Utilizes known variant sites to avoid mistaking true variants for sequencing errors.
     - Generates a recalibration table capturing these error patterns.

2. **ApplyBQSR:**
   - **Function:** Adjusts the base quality scores in the BAM file based on the recalibration table.
   - **Process:**
     - Applies the statistical models from `recal_data.table` to modify the quality scores.
     - Outputs a new BAM file (`recal_reads.bam`) with updated quality scores that more accurately reflect the true base quality.

### **Best Practices**

- **Use High-Quality Known Sites:**
  - Ensure that the `--known-sites` VCF file is well-curated and corresponds to the same reference genome build.
  
- **Sufficient Read Depth:**
  - BQSR requires adequate read depth to model error patterns effectively. Low coverage may lead to inaccurate recalibration.
  
- **Consistency:**
  - Maintain consistent reference genomes and known sites files across all steps to prevent mismatches and errors.
  
- **Resource Allocation:**
  - BQSR can be computationally intensive. Allocate sufficient memory and CPU resources, especially for large datasets.

### **Associated File Formats with Examples**

1. **Recalibration Table (`recal_data.table`):**

   **Example Content:**

   ```
   ReadGroup	QualityScore	ReportedBase	ComputedErrorRate	...
   sample_1	30	A	0.0012	...
   sample_1	30	C	0.0008	...
   ...
   ```

   - **Description:** Contains statistical models of error rates based on various covariates like quality score, base, and read group.
   - **Usage:** Used by `ApplyBQSR` to adjust base quality scores.

2. **Recalibrated BAM (`recal_reads.bam`):**

   - **Description:** A BAM file identical to the input (`dedup_reads.bam`) but with adjusted base quality scores.
   - **Usage:** Serves as the input for variant calling tools, ensuring higher accuracy in variant detection.

### **Visual Representation**

![BQSR Workflow](https://gatk.broadinstitute.org/hc/article_attachments/360052033071/2020_12_18_06_57_56_Base_Recalibration.png)

*Figure: BQSR Workflow Overview*

### **Conclusion of Step 10**

Base Quality Score Recalibration is a vital step to enhance the reliability of your variant calls. By correcting systematic errors in base quality scores, BQSR ensures that downstream analyses are based on the most accurate representation of your sequencing data.

---

## **Step 11: Variant Calling with GATK HaplotypeCaller**

### **Overview**

**Variant Calling** is the process of identifying genetic variants (such as Single Nucleotide Polymorphisms - SNPs and Insertions/Deletions - INDELs) from sequencing data. The **GATK HaplotypeCaller** is a widely used tool for this purpose, employing sophisticated algorithms to detect variants with high accuracy.

### **Why Variant Calling Is Important**

1. **Identifies Genetic Diversity:**
   - Reveals variations in the genome that may contribute to phenotypic diversity or disease.
   
2. **Enables Downstream Analyses:**
   - Provides the foundation for studies like Genome-Wide Association Studies (GWAS), population genetics, and functional genomics.
   
3. **Clinical Applications:**
   - Essential for personalized medicine, including identifying pathogenic variants for disease diagnosis and treatment.

### **Associated File Formats**

- **BAM (`.bam`):**
  - **Description:** Input file containing aligned and recalibrated sequencing reads.
  - **Usage:** Input for HaplotypeCaller.
  
- **VCF (`.vcf` or `.vcf.gz`):**
  - **Description:** Output file containing called genetic variants.
  - **Usage:** Input for variant annotation and downstream analyses.

### **Command and Parameters**

```bash
gatk HaplotypeCaller \
    -R ~/NGS_course/data/references/$ref \
    -I recal_reads.bam \
    -O ~/NGS_course/results/variants/raw_variants.vcf
```

**Breaking Down the Command:**

- **`gatk HaplotypeCaller`:**
  - Invokes the HaplotypeCaller tool from GATK.

- **Parameters:**
  - **`-R ~/NGS_course/data/references/$ref`**
    - **Description:** Specifies the reference genome file.
    - **Importance:** Essential for aligning reads and identifying variants accurately.
  
  - **`-I recal_reads.bam`**
    - **Description:** Input BAM file after BQSR.
    - **Importance:** Provides high-quality, aligned reads for accurate variant calling.
  
  - **`-O ~/NGS_course/results/variants/raw_variants.vcf`**
    - **Description:** Defines the output VCF file to store called variants.
  
- **Optional Parameters:**
  - **`--native-pair-hmm-threads 6`**
    - **Description:** Allocates 6 CPU threads for the Pair-HMM algorithm.
    - **Importance:** Speeds up the variant calling process, especially for large datasets.
  
  - **`-ERC GVCF`**
    - **Description:** Emits a genomic VCF (gVCF) which includes information about all sites, not just variant positions.
    - **Usage:** Useful for joint genotyping in multiple samples.

**Example with Optional Parameters:**

```bash
gatk HaplotypeCaller \
    -R ~/NGS_course/data/references/hg38.fa \
    -I recal_reads.bam \
    --native-pair-hmm-threads 6 \
    -O ~/NGS_course/results/variants/raw_variants.vcf
```

### **Step-by-Step Explanation**

1. **Function of HaplotypeCaller:**
   - Performs local de novo assembly of haplotypes in regions with evidence of variation.
   - Accurately identifies SNPs and INDELs by considering the sequence context and read support.

2. **Key Processes:**
   - **Local Assembly:**
     - Constructs potential haplotypes based on the reads aligned to a region.
     - Differentiates between true variants and sequencing errors.
   
   - **Genotype Likelihoods:**
     - Calculates the probability of each possible genotype (e.g., homozygous reference, heterozygous, homozygous alternate) at each variant site.
   
   - **Variant Calling:**
     - Determines the most likely genotype for each site based on the likelihoods.
     - Assigns quality scores to each variant call, indicating confidence.

3. **Resource Considerations:**
   - **CPU and Memory:** Variant calling can be resource-intensive. Allocating sufficient CPU threads (as shown with `--native-pair-hmm-threads`) can significantly speed up the process.
   - **Parallel Processing:** For multiple samples, consider running HaplotypeCaller in parallel to optimize time.

### **Associated File Formats with Examples**

1. **Variant Call Format (VCF) File (`raw_variants.vcf`):**

   **Example Entry:**

   ```
   chr1	1302015	.	C	T	602.64	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.176;DP=76;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.93;ReadPosRankSum=0.042;SOR=0.620;ANN=T|missense_variant|MODERATE|ACAP3|ENSG00000131584|transcript|ENST00000354700.10|protein_coding|5/24|c.311G>A|p.Arg104Gln|426/3793|311/2505|104/834||
   ```

   - **Fields:**
     - **CHROM:** `chr1` (Chromosome 1)
     - **POS:** `1302015` (Position on chromosome)
     - **ID:** `.` (No identifier assigned)
     - **REF:** `C` (Reference allele)
     - **ALT:** `T` (Alternate allele)
     - **QUAL:** `602.64` (Quality score)
     - **FILTER:** `PASS` (Variant passed all filters)
     - **INFO:** Contains detailed annotations and metrics (explained in previous sections)

2. **Genomic Variant Call Format (gVCF) File (`raw_variants.g.vcf`):**

   - **Description:** An extended version of VCF that includes information about all sites, not just variants.
   - **Usage:** Facilitates joint genotyping across multiple samples, enhancing variant discovery accuracy.

   **Example Entry:**

   ```
   chr1	1302015	.	C	T	602.64	PASS	AC=1;AF=0.500;AN=2;BaseQRankSum=-0.176;DP=76;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.93;ReadPosRankSum=0.042;SOR=0.620;END=1302015;SVTYPE=SNV	.	.	.	.
   ```

   - **Additional Fields:** May include `END`, `SVTYPE`, and other annotations specific to gVCF.

### **Understanding the VCF Entry in Context**

Let's revisit the VCF entry from Step 9 to see how HaplotypeCaller outputs variants.

```
chr1    1302015 .       C       T       602.64  PASS    AC=1;AF=0.500;AN=2;BaseQRankSum=-0.176;DP=76;ExcessHet=0.0000;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=7.93;ReadPosRankSum=0.042;SOR=0.620;ANN=T|missense_variant|MODERATE|ACAP3|ENSG00000131584|transcript|ENST00000354700.10|protein_coding|5/24|c.311G>A|p.Arg104Gln|426/3793|311/2505|104/834||
```

- **Interpretation:**
  - **Chromosome & Position:** The variant is located on Chromosome 1 at position 1,302,015.
  - **Alleles:** Cytosine (C) is replaced by Thymine (T).
  - **Quality & Filter:** High confidence with a QUAL score of 602.64 and passed all filters.
  - **INFO Field:** Provides comprehensive metrics and functional annotations, indicating this is a missense variant with moderate impact on the ACAP3 gene, resulting in an amino acid change from Arginine to Glutamine at position 104.

### **Best Practices**

1. **Resource Management:**
   - Allocate sufficient memory and CPU threads to handle large datasets efficiently.
   - Example: Use `--native-pair-hmm-threads` to leverage multiple CPU cores.

2. **Parameter Optimization:**
   - While default parameters are suitable for general purposes, consider tuning parameters like `--min-base-quality-score` or `--stand-call-conf` based on your data characteristics.

3. **Pipeline Consistency:**
   - Ensure that all steps use the same reference genome build to maintain consistency.
   - Example: Use `hg38.fa` consistently across alignment, BQSR, and variant calling.

4. **Quality Control:**
   - Review variant call metrics and consider additional filtering steps (e.g., Variant Quality Score Recalibration - VQSR) to enhance variant call reliability.

5. **Documentation and Version Control:**
   - Keep track of tool versions and parameters used to ensure reproducibility.
   - Example: Document the GATK version (`GATK 4.x`) and any specific flags or settings used.

### **Visual Representation**

![HaplotypeCaller Workflow](https://gatk.broadinstitute.org/hc/article_attachments/360056910932/HaplotypeCaller.png)

*Figure: HaplotypeCaller Workflow Overview*

### **Conclusion of Step 11**

Variant calling with GATK HaplotypeCaller is a sophisticated process that leverages advanced algorithms to accurately identify genetic variations within sequencing data. By following best practices and understanding each parameter and file format involved, you can ensure that your variant calls are both accurate and reliable, laying the groundwork for meaningful downstream analyses.

---

## **Summary of Steps 10 and 11**

- **Step 10 (BQSR):** Corrects systematic errors in base quality scores, producing a recalibrated BAM file (`recal_reads.bam`) that enhances the accuracy of variant calls.

- **Step 11 (HaplotypeCaller):** Identifies genetic variants from the recalibrated BAM file, outputting a VCF file (`raw_variants.vcf`) containing SNPs and INDELs for further analysis.

---

## **Key Takeaways**

1. **Importance of BQSR:**
   - Enhances the reliability of variant calls by correcting base quality score biases.
   
2. **Functionality of HaplotypeCaller:**
   - Utilizes local assembly and sophisticated algorithms to accurately detect genetic variants.
   
3. **File Formats:**
   - Understanding BAM, VCF, and related file formats is essential for navigating the NGS data analysis pipeline.
   
4. **Best Practices:**
   - Consistency in reference genomes, thorough quality control, and resource optimization are critical for successful analyses.

---


---

## **Additional Resources**

- **GATK Documentation:**
  - [BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132-BaseRecalibrator)
  - [ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132-BaseRecalibrator)
  - [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-HaplotypeCaller)
  
- **File Format Specifications:**
  - [VCF Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
  - [SAM/BAM Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
  
- **Tutorials and Workshops:**
  - [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-HaplotypeCaller)
  - [Understanding BQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132-BaseRecalibrator)

---


### Step 9: Variant Filtering

**Objective:** Apply quality filters to the called variants to improve reliability.

**Commands:**

```bash
java -jar ~/NGS_tutorial/tools/gatk.jar VariantFiltration \
-R ~/NGS_tutorial/data/references/hg38.fa \
-V raw_variants.vcf \
--filter-name "QD_filter" \
--filter-expression "QD < 2.0" \
--filter-name "FS_filter" \
--filter-expression "FS > 60.0" \
--filter-name "MQ_filter" \
--filter-expression "MQ < 40.0" \
--filter-name "SOR_filter" \
--filter-expression "SOR > 4.0" \
-O filtered_variants.vcf
```

**Explanation:**

- **`VariantFiltration`**: Applies filters to variants based on specified criteria.

**Parameters Explained:**

- **`--filter-name "QD_filter"`**: Names the filter for variants failing the following expression.
- **`--filter-expression "QD < 2.0"`**: Quality by Depth less than 2.0.
- **`FS > 60.0`**, **`MQ < 40.0`**, **`SOR > 4.0`**: Other quality metrics.

**Quality Metrics Explained:**

- **QD (Quality by Depth):** Variant confidence normalized by depth.
- **FS (Fisher Strand):** Strand bias.
- **MQ (Mapping Quality):** Root mean square of mapping quality.
- **SOR (Strand Odds Ratio):** Strand bias using a symmetric odds ratio test.

**Why This Step is Important:**

- Removes likely false positives and low-quality variants.

**Alternative Tools:**

- **VQSR (Variant Quality Score Recalibration)** in GATK (requires large datasets).
- **bcftools filter**: For filtering VCF files.

---

### Step 10: Variant Annotation with snpEff

**Objective:** Annotate variants to predict their effects on genes and proteins.

**Command:**

```bash
java -jar ~/NGS_tutorial/tools/snpEff/snpEff.jar -v GRCh38.99 \
filtered_variants.vcf > annotated_variants.vcf
```

**Explanation:**

- **`snpEff`**: Annotates variants with information like gene function and effect.
- **`-v`**: Verbose mode.
- **`GRCh38.99`**: Database version matching the reference genome.

**Parameters Explained:**

- **`-v`**: Provides detailed logging information.
- **`GRCh38.99`**: Ensures annotations match the genome build used.
- **`filtered_variants.vcf`**: Input VCF file after filtering.
- **`annotated_variants.vcf`**: Output annotated VCF file.

**Why This Step is Important:**

- Adds biological context to variants, facilitating interpretation.

**Alternative Tools:**

- **ANNOVAR**: A tool for variant annotation.
- **VEP (Variant Effect Predictor)** from Ensembl.

---

### Step 11: Extracting Specific Variants with SnpSift

**Objective:** Extract variants of interest, such as missense variants, from the annotated VCF file.

**Command:**

```bash
java -jar ~/NGS_tutorial/tools/snpEff/SnpSift.jar filter \
"(ANN[*].EFFECT has 'missense_variant')" \
annotated_variants.vcf > missense_variants.vcf
```

**Explanation:**

- **`SnpSift filter`**: Filters VCF files based on expressions.
- **`(ANN[*].EFFECT has 'missense_variant')`**: Expression to select variants with the effect 'missense_variant'.

**Parameters Explained:**

- **`ANN[*].EFFECT`**: Accesses the 'EFFECT' field in the 'ANN' (annotations) INFO field.
- **`has 'missense_variant'`**: Checks if the effect is 'missense_variant'.
- **`annotated_variants.vcf`**: Input VCF file with annotations.
- **`missense_variants.vcf`**: Output VCF file with selected variants.

**Why This Step is Important:**

- Missense variants can alter protein function and may be associated with diseases like autism.
- Focuses analysis on potentially impactful variants.

**Alternative Tools:**

- **bcftools view**: Can filter VCF files based on INFO fields.
- **VEP filter**: For filtering based on annotations.

---

## Alternative Tools and Methods

Throughout the tutorial, we have mentioned alternative tools available for similar tasks. Here's a summary:

- **Alignment:**
  - **Bowtie2**
  - **HISAT2**
- **Variant Calling:**
  - **FreeBayes**
  - **Samtools mpileup** and **bcftools**
- **Variant Annotation:**
  - **ANNOVAR**
  - **Ensembl VEP**
- **Quality Control and Trimming:**
  - **FastP**
  - **Trim Galore**
- **Variant Filtering:**
  - **bcftools filter**
  - **GATK VQSR** (for larger datasets)

**Choosing Tools:**

- The choice of tools may depend on factors like dataset size, computational resources, specific requirements of the analysis, and personal preference.

---

## Conclusion

This tutorial has guided you through a comprehensive NGS data analysis pipeline, providing detailed explanations suitable for beginners. By working through each step, you have learned:

- How to download and assess raw sequencing data.
- The importance of read alignment and how to perform it.
- Post-alignment processing steps to prepare data for variant calling.
- How to perform variant calling and filtering to obtain high-confidence variants.
- Annotating variants to understand their potential biological impact.
- Extracting variants of particular interest for further analysis.

**Next Steps:**

- **Interpretation:** Analyze the missense variants to identify those that may be associated with autism.
- **Visualization:** Use tools like IGV (Integrative Genomics Viewer) to visualize variants in the genome context.
- **Functional Analysis:** Investigate the genes affected by the variants for known associations with autism.
- **Validation:** Consider experimental validation of significant variants.

---

## Additional Resources

- **GATK Best Practices:** [https://gatk.broadinstitute.org/hc/en-us/articles/360035535912](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912)
- **snpEff Manual:** [http://snpeff.sourceforge.net/SnpEff_manual.html](http://snpeff.sourceforge.net/SnpEff_manual.html)
- **SnpSift Documentation:** [http://snpeff.sourceforge.net/SnpSift.html](http://snpeff.sourceforge.net/SnpSift.html)
- **Biostars Community:** [https://www.biostars.org/](https://www.biostars.org/) - A helpful forum for bioinformatics questions.
- **SeqAnswers Forum:** [http://seqanswers.com/](http://seqanswers.com/) - Another community resource.

---

**Note:** Bioinformatics is a rapidly evolving field. Always check for the latest versions of tools and updated best practices. The versions used in this tutorial are:

- **BWA:** Version 0.7.17 or higher
- **Samtools:** Version 1.9 or higher
- **Picard Tools:** Version 2.18.1 or higher
- **GATK:** Version 4.x
- **snpEff:** Version 5.0 or higher

**Disclaimer:** Ensure you have appropriate permissions and comply with any data usage policies when working with human genomic data.

---

**Happy Analyzing!**

