WES Exome Sequencing Analysis Pipeline (Nextflow DSL2)
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

🧬 2. Workflow Diagram
<p align="center"> <img src="workflow.svg" width="850"> </p>

3. Key Features
✔ End-to-End Automated WES Pipeline
From FASTQ → TSV with zero manual intervention.
✔ Nextflow DSL2 Modular Architecture
Each step is a standalone module (modules/*.nf), easy to edit or extend.
✔ Massive Performance Improvements
Multi-sample parallelism
Scatter execution for variant calling
Parallel annotation (1000G + ANNOVAR)
See benchmark section below.
✔ Container Support
Dockerfile included
Singularity definition included
✔ HPC/SLURM Support
Preconfigured profile in nextflow.config
Automatic container binding
✔ Reproducible & Resume-capable
-resume supported
Deterministic outputs
✔ Continuous Integration
GitHub Actions CI pipeline
Validates syntax + workflow integrity
⚙️ 4. Installation & Requirements
🔧 Software Dependencies
| Software           | Version | Purpose                         |
| ------------------ | ------- | ------------------------------- |
| **Nextflow**       | ≥ 23.x  | Workflow engine                 |
| **Java**           | ≥ 11    | Needed for GATK + Nextflow      |
| **fastp**          | Latest  | QC + Trimming                   |
| **BWA-MEM2**       | Latest  | Alignment                       |
| **Samtools**       | ≥ 1.13  | BAM processing                  |
| **GATK**           | 4.x     | Variant calling + preprocessing |
| **SnpEff/SnpSift** | Latest  | Annotation                      |
| **ANNOVAR**        | Latest  | Functional annotation           |

📁 Reference Files Required
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

