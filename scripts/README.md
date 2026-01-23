# Slurm Benchmark Script

A standalone Python script that collects Slurm job statistics for Isabl analyses and exports them to CSV for benchmarking and resource analysis.

## Features

- Collects job metrics from Slurm using `seff` and `sacct` commands
- Combines Slurm stats with Isabl analysis metadata
- Supports filtering by status, application, project, or specific analysis IDs
- Parallel processing for faster data collection
- Exports to CSV with pandas

## Requirements

- Python 3.8+
- Access to a Slurm cluster (requires `seff` and/or `sacct` commands)
- Isabl CLI configured with API access

## Setup

### 1. Create a virtual environment

```bash
# Create a new virtual environment
python3 -m venv venv

# Activate the virtual environment
# On Linux/macOS:
source venv/bin/activate

# On Windows:
# venv\Scripts\activate
```

### 2. Install dependencies

```bash
# Install isabl_cli from GitHub
pip install git+https://github.com/papaemmelab/isabl_cli.git

# Install pandas (required for data processing)
pip install pandas

# Install click (usually included with isabl_cli, but just in case)
pip install click
```

### 3. Configure Isabl CLI

Make sure your Isabl CLI is configured with the correct API URL and credentials:

```bash
# Set the API URL (replace with your Isabl instance)
export ISABL_API_URL="https://your-isabl-instance.com/api/v1/"

# If authentication is required, set your token
export ISABL_API_TOKEN="your-api-token"
```

Or configure via the Isabl settings file (`~/.isabl/settings.yaml`).

## Usage

### Basic Usage

```bash
# Collect stats for all analyses (default: hybrid mode with seff + sacct)
python slurm_benchmark.py -o benchmark.csv

# Filter by status
python slurm_benchmark.py -o benchmark.csv --status SUCCEEDED

# Filter by application
python slurm_benchmark.py -o benchmark.csv --application 123

# Filter by project
python slurm_benchmark.py -o benchmark.csv --project 456

# Specific analysis IDs
python slurm_benchmark.py -o benchmark.csv --analysis-pk 1001 --analysis-pk 1002
```

### Performance Options

```bash
# Use seff only (faster, but no timestamps/MaxRSS columns)
python slurm_benchmark.py -o benchmark.csv --seff-only

# Increase parallel workers (default: 4)
python slurm_benchmark.py -o benchmark.csv --workers 8

# Limit number of analyses (useful for testing)
python slurm_benchmark.py -o benchmark.csv --limit 10 --verbose
```

### All Options

```
Options:
  -o, --output PATH        Output CSV file path (required)
  --status TEXT            Filter by analysis status (e.g., SUCCEEDED, FAILED)
  --application INT        Filter by application PK (can specify multiple)
  --project INT            Filter by project PK (can specify multiple)
  --analysis-pk INT        Specific analysis PKs (can specify multiple)
  --seff-only              Only use seff (skip sacct for timestamps)
  -w, --workers INT        Number of parallel workers (default: 4)
  --limit INT              Limit number of analyses (for testing)
  -v, --verbose            Verbose output with summary statistics
  --help                   Show help message
```

## Output Columns

### Isabl Metadata
- `analysis_pk` - Analysis primary key
- `analysis_status` - Analysis status (SUCCEEDED, FAILED, etc.)
- `application_pk`, `application_name`, `application_version` - Application details
- `ran_by` - User who ran the analysis
- `created` - Creation timestamp
- `target_count` - Number of target experiments
- `target_system_ids` - Comma-separated target system IDs
- `project_pks` - Comma-separated project primary keys
- `storage_url` - Path to analysis results

### Slurm Job Stats
- `job_id` - Slurm job ID
- `cluster` - Slurm cluster name
- `partition` - Slurm partition (from sacct)
- `state` - Job state (COMPLETED, FAILED, etc.)
- `exit_code` - Job exit code

### Requested Resources (from sacct)
- `req_cpus` - Requested CPUs
- `req_mem`, `req_mem_bytes` - Requested memory
- `req_nodes` - Requested nodes
- `timelimit`, `timelimit_seconds` - Requested time limit

### Time Metrics
- `wall_clock_time`, `wall_clock_seconds` - Total elapsed time
- `cpu_utilized`, `cpu_utilized_seconds` - CPU time used
- `cpu_efficiency_pct` - CPU efficiency percentage
- `submit_time`, `start_time`, `end_time` - Job timestamps (from sacct)

### Allocated/Used Resources
- `cores`, `nodes` - Allocated resources
- `memory_utilized`, `memory_utilized_bytes` - Memory used (from seff)
- `memory_efficiency_pct` - Memory efficiency percentage
- `max_rss`, `max_rss_bytes` - Maximum RSS (from sacct)

### Derived Metrics
- `core_hours` - wall_clock_hours × cores
- `cpu_hours` - Actual CPU time in hours

## Data Collection Modes

### Hybrid Mode (Default)

Combines data from both `seff` and `sacct`:
- **seff**: CPU/Memory efficiency percentages (what seff does best)
- **sacct**: Partition, timestamps, MaxRSS (fields seff doesn't provide)

### Seff-Only Mode (`--seff-only`)

Uses only `seff` command. Faster but missing:
- `partition`
- `submit_time`, `start_time`, `end_time`
- `max_rss`, `max_rss_bytes`
- `req_cpus`, `req_mem`, `req_mem_bytes`, `req_nodes`
- `timelimit`, `timelimit_seconds`

## Troubleshooting

### "Slurm is not available on this system"

This script must be run on a machine with Slurm client tools installed (typically a cluster head node or login node). The `seff` and `sacct` commands must be in your PATH.

### "No analyses found"

Check your filters and ensure:
- The Isabl API is accessible
- Your authentication is configured correctly
- Analyses exist matching your filter criteria

### Empty columns

- **partition, timestamps, max_rss empty**: These come from `sacct`. Make sure you're not using `--seff-only` mode.
- **efficiency columns empty**: These come from `seff`. The job may still be running or the data may not be available.

## Example Workflow

```bash
# 1. Activate environment
source venv/bin/activate

# 2. Test with a small sample
python slurm_benchmark.py -o test.csv --limit 5 --verbose

# 3. Run full benchmark for succeeded analyses
python slurm_benchmark.py -o benchmark_succeeded.csv --status SUCCEEDED --verbose

# 4. Analyze in Python/Jupyter
python -c "
import pandas as pd
df = pd.read_csv('benchmark_succeeded.csv')
print(f'Total analyses: {df.analysis_pk.nunique()}')
print(f'Total core-hours: {df.core_hours.sum():.2f}')
print(f'Avg CPU efficiency: {df.cpu_efficiency_pct.mean():.1f}%')
print(df.groupby('application_name')['core_hours'].sum().sort_values(ascending=False))
"
```
