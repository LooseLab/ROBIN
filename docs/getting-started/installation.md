# Installation

Install ROBIN from this repository, create a fresh environment, install the package and download the required model assets.

For running ROBIN after installation, continue to the [Quickstart](quickstart.md).

## Before you start

| You need | Notes |
|----------|--------|
| **Git** | Used to clone ROBIN and its submodules. Git LFS is not required for the current repository. |
| **Conda** | Miniconda or Anaconda. |
| **Python 3.12** | Provided by the `robin` conda environment. |

**Recommended system memory:** 64 GB RAM or more for typical production-scale use. CPU/GPU requirements depend on the Oxford Nanopore sequencing setup. Docker is optional for some downstream analysis paths.

## 1. Clone the repository

```bash
git clone --recursive https://github.com/LooseLab/ROBIN.git
cd ROBIN
```

If the repository was cloned without submodules:

```bash
git submodule update --init --recursive
```

## 2. Create the conda environment

ROBIN uses the environment definition in `robin.yml`:

```bash
conda env create -f robin.yml
conda activate robin
```

Use a **fresh environment** for this codebase rather than reusing an environment from an older ROBIN release or another project.

### If the `robin` environment already exists

`robin.yml` defines `name: robin`. Choose one of the following:

| Approach | Command |
|----------|---------|
| Update in place | `conda env update -n robin -f robin.yml --prune` |
| Remove and recreate | `conda env remove -n robin` followed by `conda env create -f robin.yml` |
| Use another name | `conda env create -f robin.yml -n robin_littlejohn` |

Then activate the environment you intend to use.

### Linux `libstdc++` / `CXXABI_1.3.15` problems

If native libraries resolve against the system `libstdc++` instead of conda's copy, apply the Linux extras environment:

```bash
conda env update -n robin -f robin_linux_extras.yml
```

## 3. Install ROBIN

From the repository root:

```bash
pip install -e .
```

This installs the `robin` command-line interface from your working tree.

Optional extras are available for features that require additional packages:

| Extra | Command | Use |
|-------|---------|-----|
| GUI | `pip install -e '.[gui]'` | NiceGUI support where not already supplied by the environment |
| MinKNOW API | `pip install -e '.[minknow]'` | Programmatic MinKNOW integration |

See [MinKNOW configuration](minknow-configuration.md) for instrument integration details.

## 4. Download model and annotation assets

```bash
robin utils update-models
robin utils update-clinvar
```

Model assets are resolved from the public sources defined in ROBIN's asset manifest and are checksum-verified.

To replace existing model downloads:

```bash
robin utils update-models --overwrite
```

## 5. Verify the installation

```bash
robin --help
robin list-job-types
```

If both commands run successfully, continue to the [Quickstart](quickstart.md).

## Troubleshooting

| Problem | Check |
|---------|-------|
| Missing submodules | `git submodule update --init --recursive` |
| Model or ClinVar download failure | Check network access and retry the relevant `robin utils` command |
| Wrong environment | Use `conda env list` and activate the environment created from `robin.yml` |
| Linux native-library error | Apply `robin_linux_extras.yml` as described above |

For operational problems after ROBIN starts, use [Troubleshooting](../using-robin/troubleshooting.md).
