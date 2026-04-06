# `robin utils`

Utility subcommands for **assets**, **reference/panel staging**, and **post-hoc** MGMT summaries. The group is separate from the live [`workflow`](workflow.md) so you can update models or prepare files without starting the orchestrator.

```bash
robin utils --help
```

Subcommands documented below:

| Command | Purpose |
|---------|---------|
| [`update-models`](#robin-utils-update-models) | Download or refresh model files from the assets manifest. |
| [`update-clinvar`](#robin-utils-update-clinvar) | Update ClinVar resources used for annotation paths. |
| [`sequencing-files`](#robin-utils-sequencing-files) | Copy panel BED(s) and download or link a GRCh38 FASTA for MinKNOW / alignment. |
| [`mgmt`](#robin-utils-mgmt) | Tabulate MGMT CpG methylation from existing `mgmt_sorted.bam` outputs. |

---

## `robin utils update-models`

Downloads or refreshes **model** files (checksums via the manifest). Typical after install or when assets change.

```bash
robin utils update-models
robin utils update-models --overwrite
```

| Option | Description |
|--------|-------------|
| `--models-dir` | Override output directory (default: ROBIN models directory). |
| `--manifest` | Path to `assets.json` (default: packaged manifest or `ROBIN_ASSETS_MANIFEST`). |
| `--overwrite` | Replace existing files. |

For **private** GitHub-hosted assets, set **`GITHUB_TOKEN`** before running.

More context: [Installation — models](../getting-started/installation.md#step-4-download-models-and-clinvar).

---

## `robin utils update-clinvar`

Updates **ClinVar** in ROBIN’s resource tree (best-effort “newer NCBI drop” behaviour).

```bash
robin utils update-clinvar
```

Pair with `update-models` when setting up a new environment.

---

## `robin utils sequencing-files`

Stages **panel BED** files and a **reference FASTA** into one folder for **MinKNOW** (adaptive sampling + alignment) and for passing the same path to **`robin workflow --reference`**.

```bash
robin utils sequencing-files --panel rCNS2 --output-dir ~/references/robin_ref
```

| Option | Description |
|--------|-------------|
| `-p` / `--panel` | **Required.** Panel name (same as `--target-panel`). |
| `-r` / `--reference` | URL to download or local `.fa` / `.fa.gz`. Default: NCBI **GRCh38 no-alt** (large download). |
| `-o` / `--output-dir` | Output directory (default: `./reference_files`). |
| `-y` / `--yes` | Non-interactive confirmation. |

The command copies packaged **`{panel}_panel_name_uniq.bed`** and, when present, **`{panel}_panel_source.bed`**, then materializes the reference.

See also: [Quickstart — Reference and panel BEDs](../getting-started/quickstart.md#reference-and-panel-beds-one-command).

---

## `robin utils mgmt`

Summarises **MGMT** CpG methylation from **`mgmt_sorted.bam`** files already produced under a ROBIN output tree (offline inspection / QC).

```bash
robin utils mgmt /path/to/robin_output --recursive -o mgmt_summary.tsv
```

| Option | Description |
|--------|-------------|
| `OUTPUT_DIR` | Root directory to search (required). |
| `-r` / `--recursive` | Find `mgmt_sorted.bam` recursively. |
| `-o` / `--out` | TSV path; default `-` (stdout). |

---

## Related

- [Installation](../getting-started/installation.md)  
- [Quickstart](../getting-started/quickstart.md)  
