# Installation

This guide covers installing ROBIN and its dependencies from source for **this repository** ([ROBIN](https://github.com/LooseLab/ROBIN)).

For **running** workflows (MinKNOW settings, BAM read limits, alignment, and CLI usage), see the repository **[README](https://github.com/LooseLab/ROBIN/blob/main/README.md)**—especially *BAM read limit and MinKNOW settings* and *Usage*. This page focuses on **installation**.

---

## ROBIN with Little John

### Prerequisites

- Git and [Git LFS](https://git-lfs.com/)
- [Conda](https://docs.conda.io/) (Miniconda or Anaconda)
- Python **3.12** (supplied by the `robin` conda environment)

Recommended:

- **64 GB RAM** or more for production runs
- CPU/GPU per Oxford Nanopore guidance for your sequencer

Docker is optional (some third-party or containerized toolchains may use it).

### Step 1: Clone the repository

Clone including submodules:

```bash
git clone --recursive https://github.com/LooseLab/ROBIN.git
cd ROBIN
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

### Step 2: Create the conda environment

Use the environment file at the repository root:

| File | Purpose |
|------|---------|
| **`robin.yml`** | Main environment: Python 3.12, scientific stack, bioinformatics tools, R/Bioconductor (Linux and macOS). |
| **`robin_linux_extras.yml`** | **Linux only**, optional: merge after creating the env if you hit `libstdc++` / `CXXABI_1.3.15` issues (see README *Common issues*). |

```bash
conda env create -f robin.yml
conda activate robin
```

Notes:

- `pyproject.toml` requires Python `>=3.12`; this environment matches that.
- Prefer a **fresh** environment for this codebase rather than reusing an older ROBIN env from past releases.

#### If the `robin` environment already exists

`robin.yml` sets **`name: robin`**. If `conda env create -f robin.yml` reports the environment already exists:

**Option A — Update in place**

```bash
conda env update -n robin -f robin.yml --prune
conda activate robin
```

**Option B — Remove and recreate** (clearest if the old env is stale or mixed)

```bash
conda deactivate
conda env remove -n robin
conda env create -f robin.yml
conda activate robin
```

**Option C — New name**

```bash
conda env create -f robin.yml -n robin_littlejohn
conda activate robin_littlejohn
```

On **Linux**, if you see `CXXABI_1.3.15` / wrong `libstdc++`:

```bash
conda env update -n robin -f robin_linux_extras.yml
```

(See README *Common issues* for details.)

### Step 3: Install ROBIN (editable)

From the repository root:

```bash
pip install -e .
```

This installs the `robin` CLI from your working tree.

### Step 4: Download models and ClinVar

After installation, fetch bundled assets (SHA256-verified; set `GITHUB_TOKEN` if assets are on private GitHub):

```bash
robin utils update-models
robin utils update-clinvar
```

Private GitHub:

```bash
export GITHUB_TOKEN=your_personal_access_token
robin utils update-models
```

Force re-download models if needed:

```bash
robin utils update-models --overwrite
```

### Step 5: Verify

```bash
robin --help
robin list-job-types
```

### Troubleshooting

| Issue | What to do |
|-------|------------|
| Missing submodules | `git submodule update --init --recursive` |
| Model / ClinVar download failures | Set `GITHUB_TOKEN` if required; retry with `robin utils update-models --overwrite`; run `robin utils update-clinvar` |
| Wrong conda env | `conda env list` and activate the env you created from `robin.yml` |

### Next steps

- [Quickstart](quickstart.md)
- [README — Usage](https://github.com/LooseLab/ROBIN/blob/main/README.md#usage)
