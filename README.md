# WES Exome Sequencing Analysis Pipeline

This repository contains a Nextflow pipeline for Whole Exome Sequencing (WES) data analysis. The pipeline starts from raw FASTQ files and produces an annotated and filtered list of variants.

![Nextflow](https://img.shields.io/badge/Nextflow-v23.04.0-brightgreen)
![Workflow Status](https://github.com/imrobintomar/WES_nextflow_pipline/actions/workflows/ci.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

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



<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 700" width="1200" height="700">
<rect x="40" y="40" width="220" height="70" class="box"/>
<text x="150" y="82" text-anchor="middle" class="title">FASTQ Input</text>
<text x="150" y="100" text-anchor="middle" class="label">sample_1/2.fastq.gz</text>


<rect x="40" y="140" width="220" height="70" class="box"/>
<text x="150" y="182" text-anchor="middle" class="title">fastp</text>
<text x="150" y="200" text-anchor="middle" class="label">QC & trimmed FASTQ</text>


<rect x="40" y="240" width="220" height="70" class="box"/>
<text x="150" y="282" text-anchor="middle" class="title">BWA-MEM2</text>
<text x="150" y="300" text-anchor="middle" class="label">SAM</text>


<!-- center column -->
<rect x="420" y="40" width="320" height="70" class="box"/>
<text x="580" y="82" text-anchor="middle" class="title">Sort & MarkDuplicates</text>
<text x="580" y="100" text-anchor="middle" class="label">coordinate BAM + metrics</text>


<rect x="420" y="140" width="320" height="70" class="box"/>
<text x="580" y="182" text-anchor="middle" class="title">BQSR</text>
<text x="580" y="200" text-anchor="middle" class="label">recalibrated BAM</text>


<rect x="420" y="240" width="320" height="70" class="box"/>
<text x="580" y="282" text-anchor="middle" class="title">HaplotypeCaller (scatter)</text>
<text x="580" y="300" text-anchor="middle" class="label">per-chromosome VCFs</text>


<!-- right column -->
<rect x="820" y="140" width="300" height="70" class="box"/>
<text x="970" y="182" text-anchor="middle" class="title">Merge VCF</text>
<text x="970" y="200" text-anchor="middle" class="label">per-sample merged VCF</text>


<rect x="820" y="240" width="300" height="70" class="box"/>
<text x="970" y="282" text-anchor="middle" class="title">1000G / SnpSift</text>
<text x="970" y="300" text-anchor="middle" class="label">population freq annotation</text>


<rect x="820" y="340" width="300" height="70" class="box"/>
<text x="970" y="382" text-anchor="middle" class="title">ANNOVAR</text>
<text x="970" y="400" text-anchor="middle" class="label">multi-db annotation</text>


<rect x="820" y="440" width="300" height="70" class="box"/>
<text x="970" y="482" text-anchor="middle" class="title">Filtering & TSV</text>
<text x="970" y="500" text-anchor="middle" class="label">final clinical TSV</text>


<!-- arrows -->
<line x1="260" y1="75" x2="420" y2="75" class="arrow" />
<line x1="260" y1="175" x2="420" y2="175" class="arrow" />
<line x1="260" y1="275" x2="420" y2="275" class="arrow" />


<line x1="740" y1="275" x2="820" y2="275" class="arrow" />
<line x1="740" y1="185" x2="820" y2="185" class="arrow" />


<line x1="740" y1="375" x2="820" y2="375" class="arrow" />
<line x1="740" y1="475" x2="820" y2="475" class="arrow" />


<!-- labels -->
<text x="620" y="35" text-anchor="middle" class="label">Core preprocessing & scatter</text>


</svg>


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