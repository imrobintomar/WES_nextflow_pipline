# WES Exome Sequencing Analysis Pipeline

This repository contains a Nextflow pipeline for Whole Exome Sequencing (WES) data analysis. The pipeline starts from raw FASTQ files and produces an annotated and filtered list of variants.

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