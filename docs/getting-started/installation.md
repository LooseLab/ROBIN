# Installation

This guide will help you install ROBIN and its dependencies.

## Prerequisites

- Docker (required for running third-party tools in containers)
- Python 3.8 or higher
- Git
- Git LFS
- Conda (recommended for installation)

## Installation Methods

### Recommended Installation (Using Conda)

The recommended way to install ROBIN is using the provided conda environment files:

1. For Linux users:
```bash
conda env create -f robin.yml
```

2. For macOS users:
```bash
conda env create -f robin_osx.yml
```

3. Activate the environment:
```bash
conda activate robin
```

4. For macOS users, you will need to install the GenomicRanges R package. Launch R and run:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
```

5. Finally, install the package:
```bash
pip install .
```

### From Source

To install from source:

1. Clone the repository and initialize required components:
```bash
git clone https://github.com/LooseLab/robin/
cd ROBIN
git lfs install
git lfs pull
git submodule update --init --recursive
```

2. Create and activate the conda environment as described above.

3. Install the package:
```bash
pip install .
```

## Dependencies

ROBIN has several dependencies that will be automatically installed through the conda environment:

- R packages (including GenomicRanges for macOS users)
- Python packages
- readfish
- ont-pyguppy-client-lib (Linux only)
- And other required packages (see robin.yml or robin_osx.yml for full list)

## Verification

To verify your installation, run:

```bash
robin --version
```

## Troubleshooting

### Common Issues

1. **Docker Requirements**
   - Ensure Docker is installed and running on your system
   - Verify Docker has sufficient permissions to run containers

2. **Git LFS Issues**
   - If you encounter issues with large files, ensure Git LFS is properly installed and initialized
   - Run `git lfs install` if you haven't already

3. **Permission Errors**
   - If you encounter permission errors during installation, ensure you have appropriate permissions for Docker and the installation directory

4. **macOS Specific Issues**
   - If you encounter issues with R packages on macOS, ensure you've installed GenomicRanges as described above
   - Some tools may have different requirements on macOS vs Linux

## Next Steps

After installation, proceed to the [Quick Start Guide](quickstart.md) to begin using ROBIN. 