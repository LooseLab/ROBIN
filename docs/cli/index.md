# Command-line reference

ROBIN exposes a **`robin`** CLI (Click). This section documents each command group. For a first run, see [Installation](../getting-started/installation.md) and [Quickstart](../getting-started/quickstart.md).

## Command tree

| Command | Purpose |
|---------|---------|
| [`robin workflow`](workflow.md) | Watch a directory of BAMs and run the selected analysis pipeline (Ray, NiceGUI optional). |
| [`robin list-job-types`](jobs.md) | Print job types and queue layout. |
| [`robin list-panels`](panels.md) / [`add-panel`](panels.md#robin-add-panel) / [`remove-panel`](panels.md#robin-remove-panel) | Manage target gene panels. |
| [`robin utils …`](utils.md) | Models, ClinVar, reference/panel staging, MGMT summaries. |
| [`robin password …`](password.md) | Set the NiceGUI login password. |

Run **`robin --help`** and **`robin <command> --help`** for the exact options shipped with your install.

## Disclaimer prompt

Many commands that change state or start processing show a **disclaimer** and require you to type **`I agree`** exactly (research-use acknowledgment).

## Environment variables (selected)

| Variable | Effect |
|----------|--------|
| `GITHUB_TOKEN` | Used by `robin utils update-models` when assets are fetched from private GitHub. |
| `ROBIN_PROCESS_LARGE_BAMS` | If enabled, warns that large-BAM mode must not be used alongside live runs. |
| `LJ_BAM_THREADS` | Optional tuning for multi-threaded BAM handling (see [README — Performance](https://github.com/LooseLab/ROBIN/blob/main/README.md#performance)). |

## First-time startup experience

When you run **`robin workflow`**, ROBIN performs model checks, optional reference validation, the **research disclaimer** (`I agree`), and (if the GUI is enabled) **GUI password** prompts. See **[What happens at startup](../getting-started/startup.md)**.

## Web interface

For a **browser user guide** (navigation, each screen, reading results), see **[Using ROBIN](../using-robin/index.md)**.

## See also

- [README — Command reference](https://github.com/LooseLab/ROBIN/blob/main/README.md#command-reference) (repo mirror of common invocations)
