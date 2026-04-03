# ROBIN Quickstart

How to run ROBIN from this repository after [installation](installation.md): live BAM analysis with `robin workflow`, optional NiceGUI monitoring, and target panels.

For **disclaimer, GUI password, and startup order**, see **[What happens at startup](startup.md)**. For a **step-by-step guide to the web interface**, see **[Using ROBIN](../using-robin/index.md)**.

For deeper detail (CLI flags, MinKNOW, troubleshooting), use the **[README](https://github.com/LooseLab/ROBIN/blob/main/README.md)** and the **[command-line reference](../cli/index.md)**.

---

## What ROBIN expects

ROBIN consumes **aligned BAMs** from Oxford Nanopore sequencing (typically written in real time by MinKNOW):

| Expectation | Notes |
|-------------|--------|
| Basecalling | **HAC** is sufficient; SUP not required. |
| Methylation | **5hmC / 5mC** calling enabled in MinKNOW if your analyses need methylation. |
| Alignment | Done **in MinKNOW** — ROBIN does not realign reads. |
| BAM rollover | **Read-count–based** chunks. **Each BAM must be ≤ 50,000 reads**; we recommend **~50,000 reads per file**. Do **not** rely on **time-only** (e.g. hourly) rollover — see [README — BAM read limit](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings). |
| POD5 / FASTQ | Not required; you can turn them off if you only need BAM. |

### Memory (≤ 64 GB RAM)

On smaller machines, restart between long runs or after switching flow-cell positions. Dorado can hold GPU/host memory; restarting Dorado or the instrument after a run reduces out-of-memory risk.

---

## Reference and panel BEDs

Use **`robin utils sequencing-files`** to gather the **processed panel BED** (and optional **source BED** when packaged) plus a **GRCh38 reference FASTA** in one place for MinKNOW (adaptive sampling, alignment) and for **`robin workflow --reference`**. Run this after [installation](installation.md) when you need a consistent reference path.

```bash
robin utils sequencing-files --panel rCNS2 --output-dir ~/references/robin_ref
```

| Option | Purpose |
|--------|---------|
| `-p` / `--panel` | **Required.** Same names as `--target-panel` (built-in panels such as `rCNS2`, `AML`, `PanCan`; run `robin utils sequencing-files --help` for the list on your install). |
| `-r` / `--reference` | **Reference FASTA:** either an **HTTPS URL** to download, or a **local path** to `.fa` / `.fa.gz`. If omitted, ROBIN uses the default **NCBI GRCh38 no-alt analysis set** (UCSC-style contig names) — a **large** download; use `-r` to point at an existing file if you already have GRCh38. |
| `-o` / `--output-dir` | Output folder (default: **`./reference_files`** in the current directory). |
| `-y` / `--yes` | Skip the confirmation prompt (for scripts). |

The command copies the panel BED(s) from ROBIN resources and **materializes** the reference at the chosen location. Use the **same** reference file for **MinKNOW alignment** and for **`robin workflow --reference`**.

### Other `robin utils` commands

| Command | Purpose |
|---------|---------|
| `robin utils update-models` | Download or refresh **classification / model** assets (see [Installation](installation.md)). |
| `robin utils update-clinvar` | Download or refresh **ClinVar** under ROBIN resources for annotation paths. |
| `robin utils mgmt` | Summarise **MGMT** CpG methylation from existing `mgmt_sorted.bam` outputs (post-run inspection). |

Run **`robin utils --help`** for the full list.

---

## Run a workflow

Typical command:

```bash
robin workflow <data_folder> --work-dir <output_folder> \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center <center_id>
```

| Argument | Meaning |
|----------|---------|
| `<data_folder>` | Directory watched for incoming BAMs |
| `--work-dir` | Root directory for outputs |
| `-w` | Comma-separated analysis types |
| `--reference` | Reference FASTA (needed for most steps) |
| `--center` | Site label (e.g. `Sherwood`, `Auckland`) |
| `--target-panel` | Panel name, e.g. `rCNS2`, `PanCan` |

### Examples

```bash
# Broad workflow with a fixed panel
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center Sherwood \
  --target-panel rCNS2

# Smaller selection
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center Auckland \
  --target-panel PanCan

# More logging
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,cnv,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center New_York \
  --target-panel rCNS2 \
  --verbose \
  --log-level INFO
```

Point `--reference` at the same **GRCh38** FASTA you use for MinKNOW — for example the file produced under your [`sequencing-files`](#reference-and-panel-beds) output directory — not a different assembly or naming scheme.

---

## Commands

### `list-job-types`

```bash
robin list-job-types
```

Examples: preprocessing, bed_conversion, mgmt, cnv, target, fusion, sturgeon, nanodx, pannanodx, random_forest (see live output for your version).

### `workflow`

```bash
robin workflow /path/to/directory -w "<workflow_plan>" [OPTIONS]
```

Required:

- `-w` / `--workflow` — e.g. `mgmt,sturgeon`, or queue-qualified plans (`robin workflow --help`).
- `--center` — site ID.

Useful options (full list: `robin workflow --help`):

- `--work-dir` — output directory  
- `--reference` — reference FASTA  
- `--verbose`, `--log-level`, `--job-log-level`  
- `--no-process-existing` — only new files after start  
- `--deduplicate-jobs` — dedupe selected job types by sample  
- `--use-ray` / `--no-use-ray` — Ray execution  
- `--with-gui` / `--no-gui` — NiceGUI workflow monitor  

---

## Panel management

Built-in panels include **rCNS2**, **AML**, **PanCan**. Custom panels are registered from BED (at least four columns: chr, start, end, gene name(s)).

```bash
robin list-panels
robin add-panel /path/to/panel.bed MyCustomPanel
robin add-panel /path/to/panel.bed MyCustomPanel --validate-only
robin remove-panel MyCustomPanel
robin remove-panel MyCustomPanel --force
```

You cannot reuse reserved names: `rCNS2`, `AML`, `PanCan`.

---

## Behaviour and limitations

- **CNV** — Heuristic calls; **review visually** before any clinical interpretation.  
- **Stop with Ctrl+C** — Graceful shutdown is attempted but not guaranteed.  
- **Bugs / questions** — [Open an issue](https://github.com/LooseLab/ROBIN/issues).  

### Performance

- Batched processing on heavy paths  
- Optional `LJ_BAM_THREADS` for BAM threading where supported  
- Non-blocking GUI updates when the NiceGUI monitor is enabled  

### License and credits

Research use only; see **LICENSE** in the repository. ROBIN integrates tools such as Sturgeon, Rapid-CNS2, Readfish, cnv_from_bam, methylartist — see the repo and papers for full attribution.

---

## Next steps

- [Library preparation](library-preparation.md)  
- [MinKNOW configuration](minknow-configuration.md)  
- [README — Usage](https://github.com/LooseLab/ROBIN/blob/main/README.md#usage)  
