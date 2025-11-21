# WES Pipeline Web Portal Dockerfile
FROM python:3.11-slim

LABEL maintainer="WES Pipeline Team"
LABEL description="Web portal for WES Nextflow Pipeline"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Install Nextflow
ENV NXF_HOME=/root/.nextflow
RUN curl -fsSL get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/ \
    && chmod +x /usr/local/bin/nextflow

# Set working directory
WORKDIR /app

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p uploads results logs

# Set environment variables
ENV FLASK_APP=app.py
ENV FLASK_ENV=production
ENV PYTHONUNBUFFERED=1

# Expose port
EXPOSE 5000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:5000/health || exit 1

# Run with gunicorn
CMD ["gunicorn", "--bind", "0.0.0.0:5000", "--workers", "4", "--timeout", "300", "app:app"]
