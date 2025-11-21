/**
 * WES Pipeline Web Portal - Main JavaScript
 */

// Utility Functions
function formatBytes(bytes, decimals = 2) {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(decimals)) + ' ' + sizes[i];
}

function formatDate(dateString) {
    const date = new Date(dateString);
    return date.toLocaleString();
}

// Status Colors
const statusColors = {
    'completed': 'success',
    'failed': 'danger',
    'running': 'primary',
    'queued': 'warning',
    'uploaded': 'secondary'
};

// Toast Notifications
function showToast(message, type = 'info') {
    const toastContainer = document.getElementById('toastContainer') || createToastContainer();

    const toast = document.createElement('div');
    toast.className = `toast align-items-center text-white bg-${type} border-0`;
    toast.setAttribute('role', 'alert');
    toast.innerHTML = `
        <div class="d-flex">
            <div class="toast-body">${message}</div>
            <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast"></button>
        </div>
    `;

    toastContainer.appendChild(toast);
    const bsToast = new bootstrap.Toast(toast);
    bsToast.show();

    toast.addEventListener('hidden.bs.toast', () => toast.remove());
}

function createToastContainer() {
    const container = document.createElement('div');
    container.id = 'toastContainer';
    container.className = 'toast-container position-fixed bottom-0 end-0 p-3';
    document.body.appendChild(container);
    return container;
}

// API Helper
async function apiRequest(url, options = {}) {
    try {
        const response = await fetch(url, {
            headers: {
                'Accept': 'application/json',
                ...options.headers
            },
            ...options
        });

        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }

        return await response.json();
    } catch (error) {
        console.error('API Error:', error);
        showToast(`Error: ${error.message}`, 'danger');
        throw error;
    }
}

// Job Status Polling
class JobStatusPoller {
    constructor(jobId, updateCallback, interval = 5000) {
        this.jobId = jobId;
        this.updateCallback = updateCallback;
        this.interval = interval;
        this.timer = null;
    }

    start() {
        this.poll();
        this.timer = setInterval(() => this.poll(), this.interval);
    }

    stop() {
        if (this.timer) {
            clearInterval(this.timer);
            this.timer = null;
        }
    }

    async poll() {
        try {
            const data = await apiRequest(`/api/job/${this.jobId}/status`);
            this.updateCallback(data);

            if (data.status === 'completed' || data.status === 'failed') {
                this.stop();
            }
        } catch (error) {
            console.error('Polling error:', error);
        }
    }
}

// File Validation
function validateFiles(files) {
    const validExtensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.gz'];
    const errors = [];

    for (const file of files) {
        const fileName = file.name.toLowerCase();
        const isValid = validExtensions.some(ext => fileName.endsWith(ext));

        if (!isValid) {
            errors.push(`Invalid file type: ${file.name}`);
        }
    }

    return {
        valid: errors.length === 0,
        errors: errors
    };
}

// Progress Bar Animation
function animateProgress(element, targetValue, duration = 500) {
    const start = parseInt(element.style.width) || 0;
    const change = targetValue - start;
    const startTime = performance.now();

    function update(currentTime) {
        const elapsed = currentTime - startTime;
        const progress = Math.min(elapsed / duration, 1);

        const easeOut = 1 - Math.pow(1 - progress, 3);
        const currentValue = start + (change * easeOut);

        element.style.width = `${currentValue}%`;
        element.textContent = `${Math.round(currentValue)}%`;

        if (progress < 1) {
            requestAnimationFrame(update);
        }
    }

    requestAnimationFrame(update);
}

// Initialize on DOM Ready
document.addEventListener('DOMContentLoaded', function() {
    // Auto-dismiss alerts after 5 seconds
    const alerts = document.querySelectorAll('.alert-dismissible');
    alerts.forEach(alert => {
        setTimeout(() => {
            const bsAlert = bootstrap.Alert.getOrCreateInstance(alert);
            bsAlert.close();
        }, 5000);
    });

    // Add confirmation to delete buttons
    const deleteButtons = document.querySelectorAll('[data-confirm]');
    deleteButtons.forEach(button => {
        button.addEventListener('click', function(e) {
            if (!confirm(this.dataset.confirm)) {
                e.preventDefault();
            }
        });
    });

    // Enable tooltips
    const tooltips = document.querySelectorAll('[data-bs-toggle="tooltip"]');
    tooltips.forEach(tooltip => {
        new bootstrap.Tooltip(tooltip);
    });
});

// Export for use in templates
window.WESPortal = {
    formatBytes,
    formatDate,
    showToast,
    apiRequest,
    JobStatusPoller,
    validateFiles,
    animateProgress
};
