# Installation

This guide will help you install ROBIN and its dependencies from source.

## Prerequisites

- Docker (required for running third-party tools in containers)
- Python 3.8 or higher
- Git
- Git LFS
- Conda (recommended for environment management)

## Step 1: Clone the Repository

First, clone the ROBIN repository and initialize required components:

```bash
git clone https://github.com/LooseLab/robin.git
cd robin
git lfs install
git lfs pull
git submodule update --init --recursive
```

## Step 2: Set Up the Environment

Create and activate the conda environment:

- For Linux:
  ```bash
  conda env create -f robin.yml
  ```
- For macOS:
  ```bash
  conda env create -f robin_osx.yml
  ```

Activate the environment:
```bash
conda activate robin
```

## Step 3: Additional macOS Setup

If you are on macOS, you will need to install the GenomicRanges R package. Launch R and run:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```

## Step 4: Install ROBIN

Install the ROBIN Python package in editable mode (recommended for development and updates):

```bash
pip install -e .
```

## Dependencies

ROBIN's dependencies (Python, R, and third-party tools) are managed by the conda environment files (`robin.yml` or `robin_osx.yml`). See these files for the full list.

## Verification

To verify your installation, run:

```bash
robin --version
```

## Troubleshooting

### Common Issues

1. **Docker Requirements:**  
   Ensure Docker is installed and running on your system.  
   Verify Docker has sufficient permissions to run containers.

1. **Git LFS Issues:**  
   If you encounter issues with large files, ensure Git LFS is properly installed and initialized.  
   Run `git lfs install` if you haven't already.

1. **Permission Errors:**  
   If you encounter permission errors during installation, ensure you have appropriate permissions for Docker and the installation directory.

1. **macOS Specific Issues:**  
   If you encounter issues with R packages on macOS, ensure you've installed GenomicRanges as described above.

### MinKNOW API Version Mismatch

If you encounter an error such as:

```text
ImportError: cannot import name 'ReadEndReason' from 'minknow_api.statistics_pb2'
```

This usually means there is a mismatch between the version of the `minknow_api` Python package and the version of MinKNOW installed on your system. ROBIN requires that the `minknow_api` version matches the MinKNOW version running on your GridION or other device.

**Solution:**
Check your MinKNOW version.

Install the matching version of the `minknow_api` Python package. For example, if you are running MinKNOW 6.2.x, install the corresponding API version:

   ```bash
   pip install minknow-api==6.2.1
   ```

You can find more information and available versions at the [nanoporetech/minknow_api GitHub repository](https://github.com/nanoporetech/minknow_api).

For more details, see the discussion in [ROBIN issue #121](https://github.com/LooseLab/ROBIN/issues/121).

## Next Steps

- Proceed to the [Quick Start Guide](quickstart.md) to begin using ROBIN.
- Return to the [Home page](../index.md). 