# ![ROBIN_logo_small.png](src/robin/gui/images/ROBIN_logo_small.png) R.O.B.I.N

# ***This software is provided as is for research use only.***

**IMPORTANT: use a fresh conda environment from `robin.yml` (`conda activate robin`). Do not reuse conda environments created for older ROBIN releases or for other projects.**

**IMPORTANT (input BAMs): each file must contain 50,000 reads or fewer.** In MinKNOW, configure **read-count–based** BAM output—**we recommend one BAM every 50,000 reads**. **Do not** use **time-based** BAM rollover (MinKNOW’s typical default, e.g. hourly); it is unsupported and usually violates the read limit. Details: [BAM read limit and MinKNOW](#bam-read-limit-and-minknow-settings).

## Table of contents

- [About](#about)
- [Requirements](#requirements)
- [Installation](#installation)
  - [If the `robin` conda environment already exists](#if-the-robin-conda-environment-already-exists)
- [Common issues](#common-issues)
- [Usage](#usage)
  - [BAM read limit and MinKNOW](#bam-read-limit-and-minknow-settings)
- [Command reference](#command-reference)
- [Known issues and limitations](#known-issues-and-limitations)
- [Performance](#performance)
- [Dependencies](#dependencies)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## About

**ROBIN** (Rapid nanopOre Brain intraoperatIve classificatioN) is a bioinformatics workflow system for processing and analysing human oncology BAM data from Oxford Nanopore sequencing. It was published in *Neuro-Oncology*: [ROBIN paper](https://academic.oup.com/neuro-oncology/article/27/8/2035/8139084).

ROBIN provides automated preprocessing, multiple analysis pipelines, and real-time monitoring. It incorporates **LITTLE JOHN** (Lightweight Infrastructure for Task Tracking and Logging with Extensible Job Orchestration for High-throughput aNalysis), which handles orchestration and scaling behind the scenes.

**Capabilities** include methylation analysis, copy-number variation, fusion detection, classification workflows, a multi-threaded execution model, and a web-based GUI for monitoring, progress, and visualisation.

This repository is the canonical home for ROBIN—development, source code, and releases: [LooseLab/ROBIN](https://github.com/LooseLab/ROBIN).

## Requirements

| Resource | Notes |
|----------|--------|
| **RAM** | ≥ 64 GB recommended |
| **GPU** | As per ONT guidelines for adaptive sampling |
| **CPU** | As per ONT guidelines (more is generally better) |

With **LITTLE JOHN**, ROBIN can run two PromethION flow cells simultaneously on a Nanopore P2i (e.g. LSK114 or modified ultra-long protocol), subject to the above resources.

---

## Installation

Use **conda** so native and Python dependencies stay consistent. Create the environment from **`robin.yml`** (Python 3.12); older Python 3.9-era env files are removed.

For a step-by-step walkthrough, see [`docs/getting-started/installation.md`](docs/getting-started/installation.md).

### Prerequisites

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda

### Steps

1. **Clone with submodules**
   ```bash
   git clone --recursive https://github.com/LooseLab/ROBIN.git
   cd ROBIN
   ```

2. **Ensure submodules are current** (e.g. nanoDX, hv_rapidCNS2)
   ```bash
   git submodule update --init --recursive
   ```

3. **Create and activate the environment**
   ```bash
   conda env create -f robin.yml
   conda activate robin
   ```
   Linux and macOS share this file. On Linux, for `libstdc++` / `CXXABI_1.3.15` errors, see [Common issues](#common-issues).

4. **Install the package in editable mode**
   ```bash
   pip install -e .
   ```

5. **Download models and ClinVar data** (models are checksum-verified via the assets manifest)
   ```bash
   robin utils update-models
   robin utils update-clinvar
   ```
  
   To **re-download** models (e.g. after a failed partial run), use `robin utils update-models --overwrite`. Advanced: `python scripts/fetch_asset.py` and [`src/robin/resources/assets.json`](src/robin/resources/assets.json).

### If the `robin` conda environment already exists

`robin.yml` defines `name: robin`. **`conda env create -f robin.yml` will fail** if an environment with that name already exists.

- **Use a different env name** (keep your existing `robin` env untouched):

  ```bash
  conda env create -f robin.yml -n robin_littlejohn
  conda activate robin_littlejohn
  ```
  
- **Update the existing env** from the current file (keeps the name `robin`):

  ```bash
  conda env update -n robin -f robin.yml --prune
  conda activate robin
  ```

  `--prune` removes packages no longer listed in the YAML (where supported).

- **Remove and recreate** (closest to a clean install; preferred if the old env was from another project or an older ROBIN release):

  ```bash
  conda deactivate
  conda env remove -n robin
  conda env create -f robin.yml
  conda activate robin
  ```



After any of these, run **`pip install -e .`** again from the repository root with the env you intend to use.

More detail: [`docs/getting-started/installation.md` → If the `robin` environment already exists](docs/getting-started/installation.md#if-the-robin-environment-already-exists).

---

## Possible issue

### `libstdc++.so.6` / `CXXABI_1.3.15` (Linux: SciPy, ICU, native extensions)

The linker may use the **system** `libstdc++.so.6` (e.g. under `/lib/x86_64-linux-gnu/`) instead of conda’s (`libstdcxx-ng` in `$CONDA_PREFIX/lib`), producing errors such as:

`version 'CXXABI_1.3.15' not found (required by ... scipy ... or libicui18n ...)`

After `conda activate robin`, prefer the env libraries first:

```bash
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH}"
```

You can add that to your shell config (after conda init). 


---

## Usage

### BAM read limit and MinKNOW settings

**ROBIN requires that each BAM file contain 50,000 reads or fewer.** Larger files are outside the supported real-time workflow.

**MinKNOW configuration:** set BAM output to roll on **read count**, not on **time**. **Recommended:** **one BAM file every 50,000 reads** (smaller roll sizes, e.g. 10,000 reads per file, are also fine). Use whatever MinKNOW option splits or rotates BAMs by **number of reads** in each file.

**Do not** use MinKNOW’s **time-based** BAM settings (for example the default behaviour of writing a new BAM every fixed period such as one hour). That mode is **unsupported**: it tends to produce BAMs with far more than 50,000 reads and does not match how ROBIN expects data to arrive.

### What ROBIN expects from your sequencing setup

- BAMs from an Oxford Nanopore sequencer; **real-time HAC** basecalling (SUP not required).
- **5hmC / 5mC** methylation calling enabled in MinKNOW.
- **Real-time alignment in MinKNOW** — ROBIN does not realign reads.
- BAMs must respect the **[50,000-read limit](#bam-read-limit-and-minknow-settings)** and MinKNOW read-count output settings above.
- ROBIN does **not** consume POD5 or FASTQ; you can disable those outputs in MinKNOW if you wish.

### Memory and Dorado

On machines with **≤ 64 GB RAM**, restart the machine (or at least Dorado) before a heavy run. Dorado can retain memory indefinitely; on a P2i, after a run on position A and then B, restarting after position B (once basecalling finishes) is recommended.

### Example workflows

Primary pattern:

```bash
robin workflow <data_folder> --work-dir <output_folder> \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center <center_id>
```

| Argument | Meaning |
|----------|---------|
| `<data_folder>` | Directory containing BAM files |
| `--work-dir` | Output directory for results |
| `-w` / `--workflow` | Comma-separated job types (see [`list-job-types`](#list-job-types)) |
| `--reference` | Reference FASTA (required for many analyses) |
| `--center` | Site ID (e.g. `Sherwood`, `Auckland`, `New York`) |
| `--target-panel` | Panel for target/CNV/fusion (e.g. `rCNS2`, `AML`) |

More examples:

```bash
# Full analysis set with panel
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center Sherwood \
  --target-panel rCNS2

# Smaller workflow
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center Auckland \
  --target-panel AML

# Verbose logging
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,cnv,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center New_York \
  --target-panel rCNS2 \
  --verbose \
  --log-level INFO
```

Full CLI flags for `workflow` are listed under [Command reference](#command-reference).

---

## Command reference

### `list-job-types`

Lists job types by queue:

```bash
robin list-job-types
```

| Queue | Job types |
|-------|-----------|
| Preprocessing | `preprocessing` |
| BED conversion | `bed_conversion` |
| Analysis | `mgmt`, `cnv`, `target`, `fusion` |
| Classification | `sturgeon`, `nanodx`, `pannanodx` |
| Slow | `random_forest` |

### `workflow`

```bash
robin workflow /path/to/directory --workflow "workflow_plan" [OPTIONS]
```

**Commonly required**

- `--workflow`, `-w` — Plan such as `mgmt,sturgeon` or queue-style `preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon`
- `--center` — Center ID (e.g. `Sherwood`, `Auckland`, `New York`)

**Common options**

- `--work-dir`, `-d` — Output base directory
- `--reference`, `-r` — Reference genome (FASTA)
- `--verbose`, `-v` — Verbose output and traces
- `--no-process-existing` — Only watch for new files
- `--log-level` — `DEBUG` \| `INFO` \| `WARNING` \| `ERROR` (default: `ERROR`)
- `--job-log-level` — Per-job level, e.g. `preprocessing:DEBUG`, `mgmt:WARNING`
- `--deduplicate-jobs` — Deduplicate by sample ID for given types (e.g. `sturgeon`, `mgmt`)
- `--no-progress` — Disable file progress bars
- `--use-ray` / `--no-use-ray` — Ray distributed execution (default: on)
- `--with-gui` / `--no-gui` — NiceGUI monitor (default: on)

### Panel management

Built-in panels include `rCNS2`, `AML`. Custom panels are stored after you add them from a BED file.

**List panels**

```bash
robin list-panels
```

**Add a custom panel** (BED: ≥ 4 columns — chr, start, end, gene name(s); 4- or 6-column BED supported; multiple genes comma-separated in one region)

```bash
robin add-panel /path/to/your_panel.bed MyCustomPanel
robin add-panel /path/to/your_panel.bed MyCustomPanel --validate-only
```

Names must be non-empty and not reserved (`rCNS2`, `AML`).

**Remove a custom panel**

```bash
robin remove-panel MyCustomPanel
robin remove-panel MyCustomPanel --force
```

Built-in panels cannot be removed.

**Use in a workflow**

```bash
robin workflow /path/to/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion \
  --target-panel MyCustomPanel \
  --center Sherwood
```

---

## Known issues and limitations

**Release and scope**

1. This release is intended for testing; feedback is welcome.
2. **Real-time variant calling** is currently unavailable; it is planned to return - post run variant calling **IS** available.
3. All analyses must be interpreted by a qualified expert.

**Operation and data**

1. CNV calls use heuristics — verify by visual inspection.
2. **Ctrl+C** attempts graceful shutdown but may not always complete cleanly.
3. CSV export is in development and not yet reliable.
4. To **reanalyse** a dataset, remove the existing results under the ROBIN output folder first.
5. Other issues may exist — please [open an issue](https://github.com/LooseLab/ROBIN/issues) where possible.

---

## Performance

- **Batched processing** across analysis workflows
- **Memory-aware** behaviour for large or long runs
- **Non-blocking GUI** updates during analysis
- **Progress tracking** with live status

---

## Dependencies

Python package versions are declared in **`pyproject.toml`**. The **`robin.yml`** conda environment supplies the scientific stack, bioinformatics tools (e.g. samtools, bedtools), and R/Bioconductor packages used by the workflows.

**Notable Python libraries** include Click, Watchdog, pysam, pandas, NumPy, SciPy, ruptures, tqdm, Ray, and NiceGUI.

**External tools** include bedtools, samtools, and R (`Rscript`) for parts of the classification stack.

**Git submodules** (e.g. nanoDX, hv_rapidCNS2) must be initialised as in [Installation](#installation).

---

## License

***This software is provided "as is", and is for research use only.***

ROBIN is distributed under a **CC BY-NC 4.0** license. See the `LICENSE` file. That license does not override licenses of third-party tools bundled or invoked by ROBIN.

---

## Acknowledgments

Third-party tools and references:

- [Sturgeon](https://github.com/marcpaga/sturgeon)
- [Rapid-CNS2](https://link.springer.com/article/10.1007/s00401-022-02415-6)
- [Readfish](https://github.com/LooseLab/readfish)
- [cnv_from_bam](https://github.com/adoni5/cnv_from_bam)
- [methylartist](https://github.com/adamewing/methylartist)

Libraries include [Click](https://click.palletsprojects.com/), [Watchdog](https://python-watchdog.readthedocs.io/), [pysam](https://pysam.readthedocs.io/), [Ray](https://ray.io/), and [NiceGUI](https://nicegui.io/).

Thanks to everyone who contributed to these ecosystems, including colleagues in Nottingham and beyond. We are particularly grateful to Areeba Patel, Felix Sahm and colleagues for Rapid-CNS2. The list is non-exhaustive; the software is under active development.
