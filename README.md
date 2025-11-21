**WES Exome Sequencing Analysis Pipeline (Nextflow DSL2)**

 <p align="center"> <img src="https://img.shields.io/badge/Nextflow-v23.04.0-brightgreen"> <img src="https://img.shields.io/badge/Workflow_Status-Automated-success"> <img src="https://img.shields.io/badge/License-MIT-blue"> </p>

A high-performance, modular, and scalable pipeline for Whole Exome Sequencing (WES) analysis built with Nextflow DSL2.
From raw FASTQ files to filtered, annotated variants, the pipeline is optimized for speed, reproducibility, and HPC/cloud scalability.

Supports: Local, Docker, Singularity, and SLURM/HPC execution.

**🎯 1. Workflow Overview**

This pipeline automates the complete WES workflow using best-practice bioinformatics tools.
It ensures clean, reproducible, and efficient variant discovery and annotation.

**🔧 Major Pipeline Stages**

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
| **Final Filtering TSV**                | High-quality, rare, impactful variants → TSV                |

**🧬2. Workflow Diagram**

<p align="center"> <img src="workflow.svg" width="550"> </p>

**3. Key Features**
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

**▶️ 5. Usage**
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

---

## Web Portal

The WES Pipeline includes a user-friendly web portal for uploading exome sequencing files and running analyses through a browser interface.

### Web Portal Features

- **File Upload**: Drag-and-drop interface for uploading FASTQ files
- **Job Management**: Track, monitor, and manage multiple analysis jobs
- **Real-time Progress**: Live progress updates during pipeline execution
- **Results Download**: Download individual files or all results as ZIP
- **Log Viewer**: View pipeline logs in real-time

### Quick Start (Web Portal)

#### Option 1: Using Docker Compose (Recommended)

```bash
# Clone the repository
git clone https://github.com/imrobintomar/WES_nextflow_pipline.git
cd WES_nextflow_pipline

# Set environment variables (customize paths)
export REF_DATA_PATH=/path/to/reference/data
export SECRET_KEY=your-secure-secret-key

# Start the web portal
docker-compose up -d

# Access at http://localhost:5000
```

#### Option 2: Running Locally

```bash
# Navigate to web portal directory
cd web_portal

# Install dependencies
pip install -r requirements.txt

# Run development server
python app.py

# Or run with Gunicorn (production)
gunicorn --bind 0.0.0.0:5000 --workers 4 --timeout 300 app:app
```

### Web Portal Architecture

```
web_portal/
├── app.py                 # Flask application
├── requirements.txt       # Python dependencies
├── Dockerfile            # Container configuration
├── gunicorn.conf.py      # Production server config
├── templates/            # HTML templates
│   ├── base.html         # Base template
│   ├── index.html        # Upload page
│   ├── configure.html    # Job configuration
│   ├── job_status.html   # Progress & results
│   ├── jobs.html         # Job list
│   ├── logs.html         # Log viewer
│   └── about.html        # Pipeline info
├── static/
│   ├── css/style.css     # Styles
│   └── js/main.js        # JavaScript
├── uploads/              # Uploaded files
├── results/              # Analysis results
└── logs/                 # Pipeline logs
```

### Web Portal API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Home page with upload form |
| `/upload` | POST | Upload FASTQ files |
| `/job/<id>/configure` | GET | Configure job parameters |
| `/job/<id>/start` | POST | Start pipeline execution |
| `/job/<id>` | GET | View job status and results |
| `/api/job/<id>/status` | GET | JSON API for job status |
| `/job/<id>/logs` | GET | View pipeline logs |
| `/job/<id>/download/<file>` | GET | Download result file |
| `/job/<id>/download-all` | GET | Download all results (ZIP) |
| `/jobs` | GET | List all jobs |
| `/health` | GET | Health check endpoint |

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `SECRET_KEY` | (random) | Flask session secret key |
| `FLASK_ENV` | production | Flask environment |
| `REF_DATA_PATH` | /data/reference | Path to reference data |

### Production Deployment

For production deployment, we recommend:

1. **Use HTTPS**: Configure Nginx as a reverse proxy with SSL
2. **Set a strong SECRET_KEY**: Generate a secure random key
3. **Use persistent storage**: Mount volumes for uploads, results, and logs
4. **Configure resource limits**: Adjust Docker/Gunicorn workers based on server capacity

Example Nginx configuration is provided in the repository for production deployments.

