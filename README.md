**WES Exome Sequencing Analysis Pipeline (Nextflow DSL2)**

 <p align="center"> <img src="https://img.shields.io/badge/Nextflow-v23.04.0-brightgreen"> <img src="https://img.shields.io/badge/Workflow_Status-Automated-success"> <img src="https://img.shields.io/badge/License-MIT-blue"> </p>

A high-performance, modular, and scalable pipeline for Whole Exome Sequencing (WES) analysis built with Nextflow DSL2.
From raw FASTQ files to filtered, annotated variants, the pipeline is optimized for speed, reproducibility, and HPC/cloud scalability.

Supports: Local, Docker, Singularity, and SLURM/HPC execution.
🎯 1. Workflow Overview

This pipeline automates the complete WES workflow using best-practice bioinformatics tools.
It ensures clean, reproducible, and efficient variant discovery and annotation.

🔧 Major Pipeline Stages
Stage	Description

| Stage                              | Description                                                 |
| ---------------------------------- | ----------------------------------------------------------- |
| **fastp**                          | Quality control, adapter trimming, low-quality read removal |
| **BWA-MEM2**                       | FASTQ alignment to reference genome                         |
| **GATK SortSam & MarkDuplicates**  | Sorting + duplicate marking                                 |
| **GATK BQSR**                      | Base Quality Score Recalibration                            |
| **GATK HaplotypeCaller (scatter)** | Per-chromosome variant calling                              |
| **GATK MergeVcfs**                 | Unified per-sample VCF                                      |
| **SnpSift (1000G)**                | Population allele frequency annotation                      |
| **ANNOVAR**                        | Multi-database functional + clinical annotation             |
| **Final Filtering**                | High-quality, rare, impactful variants → TSV                |

🧬2. Workflow Diagram
<p align="center"> <img src="workflow.svg" width="850"> </p>

3. Key Features
**✔ End-to-End Automated WES Pipeline**

From FASTQ → TSV with zero manual intervention.

**✔ Nextflow DSL2 Modular Architecture**
Each step is a standalone module (modules/*.nf), easy to edit or extend.

**✔ Massive Performance Improvements**
Multi-sample parallelism
Scatter execution for variant calling
Parallel annotation (1000G + ANNOVAR)

**⚙️ 4. Installation & Requirements**

🔧 Software Dependencies


| Software           | Version | Purpose                         |
| ------------------ | ------- | ------------------------------- |
| **Nextflow**       | ≥ 23.x  | Workflow engine                 |
| **Java**           | ≥ 17    | Needed for GATK + Nextflow      |
| **fastp**          | Latest  | QC + Trimming                   |
| **BWA-MEM2**       | Latest  | Alignment                       |
| **Samtools**       | ≥ 1.13  | BAM processing                  |
| **GATK**           | 4.x     | Variant calling + preprocessing |
| **SnpEff/SnpSift** | Latest  | Annotation                      |
| **ANNOVAR**        | Latest  | Functional annotation           |

**📁 Reference Files Required**
    hg38 reference FASTA
    
    Mills & 1000G known sites (for BQSR)
    
    1000 Genomes chromosome VCFs
    
    ANNOVAR hg38 humandb

▶️ 5. Usage
Basic Execution

nextflow run main.nf \
  --input "/path/to/*.fastq.gz" \
  --ref "/path/to/hg38.fa" \
  --knownsites "/path/to/Mills_and_1000G.vcf.gz" \
  --output results/





## Workflow

The pipeline consists of the following steps:

1.  **Quality Control (fastp)**: Raw sequencing reads are processed with `fastp` to remove adapters, trim low-quality bases, and filter out short reads.

2.  **Alignment (BWA-MEM2)**: The cleaned reads are aligned to a reference genome using `bwa-mem2`.

3.  **Sorting and Duplicate Marking (GATK)**: The aligned SAM file is converted to BAM, sorted by coordinate, and duplicate reads are marked using `gatk SortSam` and `gatk MarkDuplicates`.
   

4.  **Base Quality Score Recalibration (GATK)**: Base quality scores are recalibrated using `gatk BaseRecalibrator` and `gatk ApplyBQSR` to correct for systematic errors.


5.  **Variant Calling (GATK HaplotypeCaller)**: Variants are called for each sample on a per-chromosome basis using `gatk HaplotypeCaller`.


6.  **Merge VCFs (GATK)**: The per-chromosome VCF files are merged into a single VCF file for each sample using `gatk MergeVcfs`.


7.  **Annotation (SnpSift & ANNOVAR)**:
    *   Variants are first annotated with allele frequencies from the 1000 Genomes Project using `SnpSift`.
    *   Further annotation is performed using `ANNOVAR` with multiple databases, including refGene, dbNSFP, ClinVar, gnomAD, and COSMIC.


8.  **Final Filtering**: The annotated VCF is converted to a tab-separated values (TSV) file and filtered using a custom `awk` script to select high-quality, rare, and potentially pathogenic variants.




## Dependencies

This pipeline requires the following software to be installed and available in your `PATH`:

*   [Nextflow](https://www.nextflow.io/)
*   [fastp](https://github.com/OpenGene/fastp)
*   [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
*   [GATK](https://gatk.broadinstitute.org/hc/en-us)
*   [SnpSift](https://pcingola.github.io/SnpEff/ss_introduction/) (part of SnpEff)
*   [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)

## Configuration

The pipeline is configured through the `nextflow.config` file. The main parameters to set are:

*   `params.input`: Path to the input FASTQ files (e.g., `"/path/to/fastqs/*.fastq.gz"`).
*   `params.ref`: Path to the reference genome FASTA file.
*   `params.knownsites`: Path to a VCF file of known variant sites for BQSR (e.g., Mills and 1000G gold standard indels).
*   `params.output`: The directory where the results will be saved.

You can also adjust the `process` selectors to configure CPU, memory, and time resources for the pipeline tasks.

## Usage

To run the pipeline, use the following command:


```bash
nextflow run main.nf -profile standard
```

Make sure you have configured the paths in `nextflow.config` before running the pipeline. The final filtered variant list will be available in the output directory specified by `params.output`.

