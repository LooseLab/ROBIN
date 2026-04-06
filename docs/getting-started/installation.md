# Installation

!!! abstract "What this page covers"
    Install ROBIN from this repository: **clone** → **conda environment** → **editable pip install** → **models & ClinVar** → **verify**.  
    For **running** workflows (MinKNOW, BAM limits, CLI examples), see the repo **[README](https://github.com/LooseLab/ROBIN/blob/main/README.md)**. This page is **install only**.

---

## Before you start

| You need | Notes |
|----------|--------|
| **Git** + **[Git LFS](https://git-lfs.com/)** | Clone includes submodules / LFS assets where used. |
| **[Conda](https://docs.conda.io/)** | Miniconda or Anaconda. |
| **Python 3.12** | Provided by the `robin` conda env (`robin.yml`). |

**Recommended for production:** **64 GB RAM** or more; CPU/GPU per Oxford Nanopore guidance for your sequencer. Docker is optional.

---

## ROBIN with Little John

### Step 1: Clone the repository

Include submodules:

```bash
git clone --recursive https://github.com/LooseLab/ROBIN.git
cd ROBIN
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

### Step 2: Create the conda environment

| File | Role |
|------|------|
| **`robin.yml`** | Main env: Python 3.12, scientific stack, bioinformatics, R/Bioconductor (Linux and macOS). |
| **`robin_linux_extras.yml`** | **Linux only**, optional: use if you hit `libstdc++` / `CXXABI_1.3.15` (see README *Common issues*). |

```bash
conda env create -f robin.yml
conda activate robin
```

- `pyproject.toml` requires Python **≥ 3.12**; this env matches that.  
- Prefer a **fresh** env for this codebase—not an older ROBIN env from past releases.

#### If the `robin` environment already exists

`robin.yml` sets **`name: robin`**. If create fails because the env exists:

| Approach | Command |
|----------|---------|
| **Update in place** | `conda env update -n robin -f robin.yml --prune` then `conda activate robin` |
| **Remove and recreate** | `conda deactivate` → `conda env remove -n robin` → `conda env create -f robin.yml` |
| **New name** | `conda env create -f robin.yml -n robin_littlejohn` → `conda activate robin_littlejohn` |

**Linux** — if you see `CXXABI_1.3.15` / wrong `libstdc++`:

```bash
conda env update -n robin -f robin_linux_extras.yml
```

(See README *Common issues*.)

### Step 3: Install ROBIN (editable)

From the repository root:

```bash
pip install -e .
```

This installs the `robin` CLI from your working tree.

### Step 4: Download models and ClinVar

Assets are SHA256-verified. Set **`GITHUB_TOKEN`** if downloads use private GitHub.

```bash
robin utils update-models
robin utils update-clinvar
```

Private GitHub:

```bash
export GITHUB_TOKEN=your_personal_access_token
robin utils update-models
```

Force re-download models:

```bash
robin utils update-models --overwrite
```

### Step 5: Verify

```bash
robin --help
robin list-job-types
```

---

## Troubleshooting (install)

| Issue | What to do |
|-------|------------|
| Missing submodules | `git submodule update --init --recursive` |
| Model / ClinVar download failures | Set `GITHUB_TOKEN` if required; retry `robin utils update-models --overwrite` and `robin utils update-clinvar` |
| Wrong conda env | `conda env list` — activate the env created from `robin.yml` |

---

## Next steps

- [Quickstart](quickstart.md) — run a workflow  
- [README — Usage](https://github.com/LooseLab/ROBIN/blob/main/README.md#usage) — deep operational detail  
