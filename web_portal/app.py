#!/usr/bin/env python3
"""
WES Nextflow Pipeline Web Portal
A Flask-based web application for whole exome sequencing analysis.
"""

import os
import uuid
import json
import subprocess
import threading
import time
from datetime import datetime
from pathlib import Path
from functools import wraps

from flask import (
    Flask, render_template, request, jsonify, send_file,
    redirect, url_for, flash, session
)
from werkzeug.utils import secure_filename

# Configuration
app = Flask(__name__)
app.secret_key = os.environ.get('SECRET_KEY', 'wes-pipeline-secret-key-change-in-production')

# Directories
BASE_DIR = Path(__file__).parent
UPLOAD_DIR = BASE_DIR / 'uploads'
RESULTS_DIR = BASE_DIR / 'results'
LOGS_DIR = BASE_DIR / 'logs'
PIPELINE_DIR = BASE_DIR.parent

# Ensure directories exist
for d in [UPLOAD_DIR, RESULTS_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Configuration
ALLOWED_EXTENSIONS = {'fastq', 'fq', 'gz', 'fastq.gz', 'fq.gz'}
MAX_CONTENT_LENGTH = 50 * 1024 * 1024 * 1024  # 50GB max upload
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH

# Job storage (in production, use a database)
JOBS_FILE = BASE_DIR / 'jobs.json'

def load_jobs():
    """Load jobs from JSON file."""
    if JOBS_FILE.exists():
        with open(JOBS_FILE, 'r') as f:
            return json.load(f)
    return {}

def save_jobs(jobs):
    """Save jobs to JSON file."""
    with open(JOBS_FILE, 'w') as f:
        json.dump(jobs, f, indent=2, default=str)

def allowed_file(filename):
    """Check if file extension is allowed."""
    if '.' not in filename:
        return False
    # Handle double extensions like .fastq.gz
    if filename.endswith('.fastq.gz') or filename.endswith('.fq.gz'):
        return True
    ext = filename.rsplit('.', 1)[1].lower()
    return ext in ALLOWED_EXTENSIONS

def get_job_status(job_id):
    """Get the current status of a job."""
    jobs = load_jobs()
    if job_id not in jobs:
        return None
    return jobs[job_id]

def update_job_status(job_id, status, progress=None, message=None, results=None):
    """Update job status."""
    jobs = load_jobs()
    if job_id in jobs:
        jobs[job_id]['status'] = status
        jobs[job_id]['updated_at'] = datetime.now().isoformat()
        if progress is not None:
            jobs[job_id]['progress'] = progress
        if message:
            jobs[job_id]['message'] = message
        if results:
            jobs[job_id]['results'] = results
        save_jobs(jobs)

def run_pipeline(job_id, input_files, config):
    """Run the Nextflow pipeline in a background thread."""
    job_dir = RESULTS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    log_file = LOGS_DIR / f"{job_id}.log"

    try:
        update_job_status(job_id, 'running', progress=5, message='Starting pipeline...')

        # Build input pattern
        upload_dir = UPLOAD_DIR / job_id
        input_pattern = str(upload_dir / "*.fastq.gz")

        # Build Nextflow command
        cmd = [
            'nextflow', 'run', str(PIPELINE_DIR / 'main.nf'),
            '--input', input_pattern,
            '--output', str(job_dir),
            '-with-report', str(job_dir / 'report.html'),
            '-with-timeline', str(job_dir / 'timeline.html'),
            '-with-dag', str(job_dir / 'dag.png'),
            '-resume'
        ]

        # Add optional parameters
        if config.get('ref'):
            cmd.extend(['--ref', config['ref']])
        if config.get('knownsites'):
            cmd.extend(['--knownsites', config['knownsites']])

        # For stub run (testing without actual processing)
        if config.get('stub_run', False):
            cmd.append('-stub-run')

        update_job_status(job_id, 'running', progress=10, message='Pipeline initialized, processing samples...')

        # Run the pipeline
        with open(log_file, 'w') as log:
            process = subprocess.Popen(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                cwd=str(PIPELINE_DIR)
            )

            # Monitor progress by checking log file
            stages = [
                ('FASTP', 15, 'Running quality control and trimming...'),
                ('BWA', 25, 'Aligning reads to reference genome...'),
                ('SORT_MARKDUP', 40, 'Sorting and marking duplicates...'),
                ('BQSR', 50, 'Applying base quality score recalibration...'),
                ('HC', 60, 'Calling variants per chromosome...'),
                ('MERGE_VCF', 75, 'Merging VCF files...'),
                ('SNP1000G', 80, 'Annotating with 1000 Genomes...'),
                ('ANNOVAR', 85, 'Running ANNOVAR annotation...'),
                ('FINAL_FILTER', 95, 'Applying final filters...')
            ]

            current_stage = 0
            while process.poll() is None:
                time.sleep(5)

                # Check log for stage progress
                if log_file.exists():
                    with open(log_file, 'r') as lf:
                        log_content = lf.read()
                        for i, (stage, prog, msg) in enumerate(stages):
                            if stage in log_content and i >= current_stage:
                                current_stage = i
                                update_job_status(job_id, 'running', progress=prog, message=msg)

            return_code = process.returncode

        if return_code == 0:
            # Find result files
            results = []
            for f in job_dir.glob('**/*'):
                if f.is_file():
                    results.append({
                        'name': f.name,
                        'path': str(f.relative_to(job_dir)),
                        'size': f.stat().st_size
                    })

            update_job_status(
                job_id, 'completed', progress=100,
                message='Pipeline completed successfully!',
                results=results
            )
        else:
            update_job_status(
                job_id, 'failed', progress=0,
                message=f'Pipeline failed with exit code {return_code}. Check logs for details.'
            )

    except Exception as e:
        update_job_status(
            job_id, 'failed', progress=0,
            message=f'Pipeline error: {str(e)}'
        )

# Routes
@app.route('/')
def index():
    """Home page with upload form."""
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    """Handle file uploads and create a new job."""
    if 'files' not in request.files:
        flash('No files selected', 'error')
        return redirect(url_for('index'))

    files = request.files.getlist('files')

    if not files or files[0].filename == '':
        flash('No files selected', 'error')
        return redirect(url_for('index'))

    # Validate files
    valid_files = []
    for file in files:
        if file and allowed_file(file.filename):
            valid_files.append(file)
        else:
            flash(f'Invalid file type: {file.filename}', 'warning')

    if not valid_files:
        flash('No valid FASTQ files found', 'error')
        return redirect(url_for('index'))

    # Create job
    job_id = str(uuid.uuid4())[:8]
    job_upload_dir = UPLOAD_DIR / job_id
    job_upload_dir.mkdir(parents=True, exist_ok=True)

    # Save files
    uploaded_files = []
    for file in valid_files:
        filename = secure_filename(file.filename)
        filepath = job_upload_dir / filename
        file.save(str(filepath))
        uploaded_files.append({
            'name': filename,
            'size': filepath.stat().st_size
        })

    # Create job entry
    jobs = load_jobs()
    jobs[job_id] = {
        'id': job_id,
        'status': 'uploaded',
        'progress': 0,
        'message': 'Files uploaded, ready to start analysis',
        'files': uploaded_files,
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat(),
        'config': {},
        'results': []
    }
    save_jobs(jobs)

    flash(f'Successfully uploaded {len(uploaded_files)} files', 'success')
    return redirect(url_for('configure_job', job_id=job_id))

@app.route('/job/<job_id>/configure')
def configure_job(job_id):
    """Configure job parameters before running."""
    job = get_job_status(job_id)
    if not job:
        flash('Job not found', 'error')
        return redirect(url_for('index'))
    return render_template('configure.html', job=job)

@app.route('/job/<job_id>/start', methods=['POST'])
def start_job(job_id):
    """Start the pipeline for a job."""
    job = get_job_status(job_id)
    if not job:
        return jsonify({'error': 'Job not found'}), 404

    if job['status'] not in ['uploaded', 'failed']:
        return jsonify({'error': 'Job cannot be started in current state'}), 400

    # Get configuration from form
    config = {
        'ref': request.form.get('ref_genome', ''),
        'knownsites': request.form.get('knownsites', ''),
        'stub_run': request.form.get('stub_run') == 'on'
    }

    # Update job config
    jobs = load_jobs()
    jobs[job_id]['config'] = config
    jobs[job_id]['status'] = 'queued'
    jobs[job_id]['message'] = 'Job queued, starting soon...'
    save_jobs(jobs)

    # Start pipeline in background thread
    thread = threading.Thread(
        target=run_pipeline,
        args=(job_id, job['files'], config)
    )
    thread.daemon = True
    thread.start()

    return redirect(url_for('job_status', job_id=job_id))

@app.route('/job/<job_id>')
def job_status(job_id):
    """View job status and results."""
    job = get_job_status(job_id)
    if not job:
        flash('Job not found', 'error')
        return redirect(url_for('index'))
    return render_template('job_status.html', job=job)

@app.route('/api/job/<job_id>/status')
def api_job_status(job_id):
    """API endpoint for job status (for polling)."""
    job = get_job_status(job_id)
    if not job:
        return jsonify({'error': 'Job not found'}), 404
    return jsonify(job)

@app.route('/job/<job_id>/logs')
def job_logs(job_id):
    """View job logs."""
    log_file = LOGS_DIR / f"{job_id}.log"
    if not log_file.exists():
        return jsonify({'logs': 'No logs available yet.'})

    with open(log_file, 'r') as f:
        logs = f.read()

    if request.headers.get('Accept') == 'application/json':
        return jsonify({'logs': logs})

    return render_template('logs.html', job_id=job_id, logs=logs)

@app.route('/job/<job_id>/download/<path:filepath>')
def download_result(job_id, filepath):
    """Download a result file."""
    job = get_job_status(job_id)
    if not job:
        flash('Job not found', 'error')
        return redirect(url_for('index'))

    file_path = RESULTS_DIR / job_id / filepath
    if not file_path.exists():
        flash('File not found', 'error')
        return redirect(url_for('job_status', job_id=job_id))

    return send_file(str(file_path), as_attachment=True)

@app.route('/job/<job_id>/download-all')
def download_all_results(job_id):
    """Download all results as a ZIP file."""
    import zipfile
    import io

    job = get_job_status(job_id)
    if not job:
        flash('Job not found', 'error')
        return redirect(url_for('index'))

    job_dir = RESULTS_DIR / job_id
    if not job_dir.exists():
        flash('Results not found', 'error')
        return redirect(url_for('job_status', job_id=job_id))

    # Create ZIP in memory
    memory_file = io.BytesIO()
    with zipfile.ZipFile(memory_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for file_path in job_dir.rglob('*'):
            if file_path.is_file():
                arcname = file_path.relative_to(job_dir)
                zf.write(file_path, arcname)

    memory_file.seek(0)
    return send_file(
        memory_file,
        mimetype='application/zip',
        as_attachment=True,
        download_name=f'wes_results_{job_id}.zip'
    )

@app.route('/jobs')
def list_jobs():
    """List all jobs."""
    jobs = load_jobs()
    # Sort by created_at descending
    sorted_jobs = sorted(
        jobs.values(),
        key=lambda x: x.get('created_at', ''),
        reverse=True
    )
    return render_template('jobs.html', jobs=sorted_jobs)

@app.route('/job/<job_id>/delete', methods=['POST'])
def delete_job(job_id):
    """Delete a job and its files."""
    import shutil

    jobs = load_jobs()
    if job_id not in jobs:
        flash('Job not found', 'error')
        return redirect(url_for('list_jobs'))

    # Don't delete running jobs
    if jobs[job_id]['status'] == 'running':
        flash('Cannot delete a running job', 'error')
        return redirect(url_for('list_jobs'))

    # Delete files
    upload_dir = UPLOAD_DIR / job_id
    results_dir = RESULTS_DIR / job_id
    log_file = LOGS_DIR / f"{job_id}.log"

    if upload_dir.exists():
        shutil.rmtree(upload_dir)
    if results_dir.exists():
        shutil.rmtree(results_dir)
    if log_file.exists():
        log_file.unlink()

    # Remove from jobs
    del jobs[job_id]
    save_jobs(jobs)

    flash('Job deleted successfully', 'success')
    return redirect(url_for('list_jobs'))

@app.route('/about')
def about():
    """About page with pipeline information."""
    return render_template('about.html')

@app.route('/health')
def health():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.now().isoformat(),
        'version': '1.0.0'
    })

# Error handlers
@app.errorhandler(413)
def too_large(e):
    flash('File too large. Maximum size is 50GB.', 'error')
    return redirect(url_for('index'))

@app.errorhandler(404)
def not_found(e):
    return render_template('404.html'), 404

@app.errorhandler(500)
def server_error(e):
    return render_template('500.html'), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
