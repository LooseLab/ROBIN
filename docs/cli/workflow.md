# `robin workflow`

This command runs the ROBIN workflow engine on BAM files under a watched directory. For a first run, start with the [Quickstart](../getting-started/quickstart.md); this page is the detailed command reference.

ROBIN can run preprocessing, BED conversion, MGMT, CNV, target analysis, fusion/ITD detection and multiple classification jobs, with optional NiceGUI monitoring.

## Synopsis

```bash
robin workflow <PATH> -w <WORKFLOW> --center <ID> --target-panel <PANEL> [OPTIONS]
```

Or load repeatable settings from TOML:

```bash
robin workflow --toml my_settings.toml
```

CLI flags override values from the TOML file when passed explicitly.

| Argument / option | Required | Description |
|-------------------|----------|-------------|
| `PATH` | Yes* | Directory containing or receiving BAM files. |
| `-t` / `--toml` | No | TOML file containing workflow settings. |
| `-w` / `--workflow` | Yes* | Comma-separated job types or legacy `queue:job` steps. |
| `--center` | Yes* | Site or study label used in outputs and reports. |
| `--target-panel` | Yes* | Built-in or custom panel name. |
| `-d` / `--work-dir` | No | Base directory for run outputs. |
| `-r` / `--reference` | No* | Reference FASTA; required by analyses that need a reference. |

\* Required on the command line or supplied through TOML where applicable.

## Configuration file

An example ships at `examples/workflow.example.toml`.

```toml
path = "empty_folder"
workflow = "cnv,fusion,target,mgmt,sturgeon,nanodx,pannanodx,random_forest"
center = "NUH"
target_panel = "rCNS2"
work_dir = "../../REF_SAMPLES"
reference = "~/references/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
log_level = "INFO"
analysis_workers = 2
with_gui = true
deduplicate_jobs = ["sturgeon", "mgmt"]
```

TOML keys use the same names as long-form CLI flags. List-valued options can be TOML arrays. `workflow` can be a comma-separated string or an array of job types.

## Workflow strings

### Simplified format — recommended

```bash
-w mgmt,sturgeon
-w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest
```

ROBIN inserts required upstream stages such as `preprocessing` and, for classifiers that require it, `bed_conversion`.

See [Job types](jobs.md) for the currently registered jobs and their queue mapping.

### Legacy explicit-queue format

Legacy workflows can specify `queue:job` pairs, for example:

```text
preprocessing:preprocessing,bed_conversion:bed_conversion,mgmt:mgmt,classification:sturgeon
```

Use the simplified format for new configurations unless explicit queue control is required.

## Sequencing input requirements

The supported real-time workflow expects BAMs that:

- are aligned before ROBIN receives them;
- use a reference compatible with the FASTA supplied through `--reference`;
- contain the required modified-base tags for methylation analyses;
- are rotated by read count rather than elapsed time;
- contain **50,000 reads or fewer per BAM**.

See [MinKNOW configuration](../getting-started/minknow-configuration.md) for instrument-specific setup.

## Execution engine

| Mode | Flag | Notes |
|------|------|-------|
| Ray | `--use-ray` | Default distributed/task-oriented execution path. |
| Threading | `--no-use-ray` | Threaded fallback execution. |

Presets such as `p2i`, `standard` and `high` adjust worker grouping and resource allocation for different hardware profiles.

## Web interface

| Option | Default | Description |
|--------|---------|-------------|
| `--with-gui` / `--no-gui` | On | Enable or disable NiceGUI monitoring. |
| `--gui-host` | `0.0.0.0` | Bind address. |
| `--gui-port` | `8081` | Port. |

The exact startup sequence, consent and initial-user behaviour are documented under [Starting ROBIN](../getting-started/startup.md).

## Logging and progress

| Option | Description |
|--------|-------------|
| `--log-level` | Global logging level. |
| `--job-log-level` | Per-job logging level, e.g. `preprocessing:DEBUG`. |
| `--verbose` / `-v` | Verbose CLI output and traces. |
| `--no-progress` | Disable progress bars. |

## File handling

| Option | Description |
|--------|-------------|
| `--no-process-existing` | Ignore BAMs already present at startup and process only newly arriving files. |
| `--no-watch` | Do not continue watching for new BAM files. |

## Ray and worker tuning

Advanced options include `--ray-num-cpus`, `--queue-priority`, `--analysis-workers`, `--preprocessing-workers`, `--bed-workers`, `--show-priorities` and `--legacy-analysis-queue`.

Use:

```bash
robin workflow --help
```

for the definitive option list in the installed version.

## Deduplication and custom commands

`--deduplicate-jobs` can restrict selected job types to a single execution per sample where appropriate.

`--commands` / `-c` can attach custom shell commands to job types. This is an advanced integration feature; verify behaviour carefully before using it in automated runs.

## Stopping a workflow

Use **Ctrl+C** to request graceful shutdown. ROBIN attempts to stop watchers, workers, Ray and the GUI cleanly, although complex runs may take a short period to terminate fully.

## Related

- [Quickstart](../getting-started/quickstart.md)
- [Starting ROBIN](../getting-started/startup.md)
- [Job types](jobs.md)
- [Panel commands](panels.md)
- [Utilities](utils.md)
