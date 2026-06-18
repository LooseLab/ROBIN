# `robin workflow`

!!! abstract "What this page covers"
    Flags and behaviour for **`robin workflow`**: synopsis, workflow string formats, Ray vs threading, GUI, logging, file handling, and links to **job types** and **quickstart**. Startup order (disclaimer, password, when the GUI appears): **[What happens at startup](../getting-started/startup.md)**.

Run the **Little John** orchestrated pipeline on BAM files under a watched directory: preprocessing, optional BED conversion, analyses (MGMT, CNV, target, fusion), classifiers (Sturgeon, NanoDX, PanNanoDX, random forest), and optional **NiceGUI** monitoring.

!!! info "Startup sequence"
    When you run this command, ROBIN checks models, optionally validates `--reference`, asks you to type **`I agree`** to the research disclaimer, then (with the GUI enabled) may prompt for the **GUI password**. On the default **Ray** path, the browser UI is only started if **`--work-dir`** is set. Details: **[What happens at startup](../getting-started/startup.md)**.

## Synopsis

```bash
robin workflow <PATH> -w <WORKFLOW> --center <ID> --target-panel <PANEL> [OPTIONS]
```

Or load settings from a TOML file:

```bash
robin workflow --toml my_settings.toml
```

CLI flags override values from the TOML file when you pass them explicitly. You can also mix both — for example, keep a site-specific config file and override the input directory per run:

```bash
robin workflow /path/to/bams --toml my_settings.toml
```

| Argument / option | Required | Description |
|-------------------|----------|-------------|
| `PATH` | Yes* | Directory containing (or receiving) BAM files. Must exist. |
| `-t` / `--toml` | No | TOML file with workflow settings (see [Configuration file](#configuration-file)). |
| `-w` / `--workflow` | Yes* | Comma-separated job types or legacy `queue:job` steps (see [Job types](jobs.md)). |
| `--center` | Yes* | Site or study label (e.g. `Sherwood`, `Auckland`) — used in outputs and reports. |
| `--target-panel` | Yes* | Panel name: built-in (`rCNS2`, `AML`, …) or custom from `robin add-panel`. Run `robin utils sequencing-files --help` to see choices on your install. |
| `-d` / `--work-dir` | No | Base directory for all run outputs. |
| `-r` / `--reference` | No* | Path to reference **FASTA**. If provided, ROBIN validates the file and ensures an index (e.g. `.fai`). Required for analyses that need a reference. |

\* Required on the command line **or** in the TOML file when using `--toml`.

## Configuration file

Use a TOML file to store repeatable run settings. An example ships with the repository at [`examples/workflow.example.toml`](https://github.com/LooseLab/ROBIN/blob/main/examples/workflow.example.toml).

```toml
path = "empty_folder"
workflow = "cnv,fusion,target,mgmt,sturgeon,nanodx,pannanodx,random_forest"
center = "NUH"
target_panel = "rCNS2"

work_dir = "../../REF_SAMPLES"
reference = "~/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

# Optional
log_level = "INFO"
analysis_workers = 2
with_gui = true
deduplicate_jobs = ["sturgeon", "mgmt"]
```

Run:

```bash
robin workflow --toml my_settings.toml
```

TOML keys use the same names as long-form CLI flags (`target_panel`, `work_dir`, `no_process_existing`, etc.). List-valued options (`deduplicate_jobs`, `job_log_level`, `commands`, `queue_priority`) can be TOML arrays. `workflow` can be a comma-separated string or a TOML array of job types.

## Workflow string formats

### Simplified (recommended)

Comma-separated **job types** only — queues are assigned automatically:

```bash
-w mgmt,sturgeon
-w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest
```

Behaviour:

- **`preprocessing`** is always inserted as the first step if missing.
- **`bed_conversion`** is inserted **before** classification jobs when needed (`sturgeon`, `nanodx`, `pannanodx`, `random_forest`).

### Legacy (explicit queues)

Use `queue:job` steps, e.g.:

```text
preprocessing:preprocessing,bed_conversion:bed_conversion,mgmt:mgmt,classification:sturgeon
```

## Execution engine

| Mode | Flag | Notes |
|------|------|--------|
| **Ray** (default) | `--use-ray` | Distributed task execution; optional Ray dashboard (`--ray-dashboard` / `--no-ray-dashboard`). |
| Threading | `--no-use-ray` | Falls back to threaded workers; analysis worker counts apply per queue mode. |

**Presets** (`--preset`): `p2i`, `standard` (default), `high` — adjust Ray actor grouping and CPU caps for different hardware (e.g. P2i vs server).

## GUI (NiceGUI)

| Option | Default | Description |
|--------|---------|-------------|
| `--with-gui` / `--no-gui` | GUI on | Launch workflow monitor in the browser. |
| `--gui-host` | `0.0.0.0` | Bind address. |
| `--gui-port` | `8081` | Port for the GUI. |

## Logging and progress

| Option | Description |
|--------|-------------|
| `--log-level` | Global level: `DEBUG`, `INFO`, `WARNING`, `ERROR` (default `ERROR`). |
| `--job-log-level` | Repeatable, e.g. `preprocessing:DEBUG`, `mgmt:WARNING`. |
| `--verbose` / `-v` | Verbose CLI and traces. |
| `--no-progress` | Disable file progress bars. |

## File handling

| Option | Description |
|--------|-------------|
| `--no-process-existing` | Only process BAMs that appear **after** startup (skip files already on disk). |
| `--no-watch` | Do not watch the directory for new files (advanced; default is to watch). |

## Ray tuning

| Option | Description |
|--------|-------------|
| `--ray-num-cpus` | Cap CPUs for Ray (auto if omitted). |
| `--queue-priority` | Repeatable, e.g. `preprocessing:10`, `mgmt:5` (Ray mode). |
| `--show-priorities` | Print queue priorities and exit. |
| `--analysis-workers` | Workers per analysis queue (default from code constant, often `1`). |
| `--preprocessing-workers`, `--bed-workers` | Workers for preprocessing and `bed_conversion`. |
| `--legacy-analysis-queue` | Single shared analysis queue instead of per-type queues (threading path). |

## Deduplication

`--deduplicate-jobs` (repeatable) — job types that should run **at most once per sample** even when multiple triggers fire (e.g. `sturgeon`, `mgmt`).

## Custom per-job shell commands

`--commands` / `-c` (repeatable) — optional mappings `job_type:shell_command` for custom steps (advanced; see `robin workflow --help`).

## Inputs and sequencing assumptions

ROBIN expects **real-time** BAMs with **≤ 50,000 reads per file** and **read-count–based** rollover in MinKNOW. See the [README](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings).

## Exit

Stop with **Ctrl+C** — the runner attempts graceful shutdown; complex runs may not always exit instantly.

## Related

- [Job types and queues](jobs.md)  
- [Panel commands](panels.md)  
- [Utilities — reference files](utils.md#robin-utils-sequencing-files)  
- [Quickstart](../getting-started/quickstart.md)  
