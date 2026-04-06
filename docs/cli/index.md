# Command-line reference

!!! abstract "What this section covers"
    The **`robin`** CLI (Click): **`workflow`**, **job types**, **panels**, **utilities**, and **GUI password**. For a first run, see [Installation](../getting-started/installation.md) and [Quickstart](../getting-started/quickstart.md).

Run **`robin --help`** and **`robin <command> --help`** for options in your build.

---

## Commands

<div class="grid" markdown>

<div class="robin-feature-card" markdown>
### [`robin workflow`](workflow.md)
Watch a BAM directory and run the selected pipeline (Ray, optional NiceGUI).
</div>

<div class="robin-feature-card" markdown>
### [`robin list-job-types`](jobs.md)
Print job types and queue layout for `-w` strings.
</div>

<div class="robin-feature-card" markdown>
### [Panels](panels.md)
`list-panels`, `add-panel`, `remove-panel` — built-in and custom BED panels.
</div>

<div class="robin-feature-card" markdown>
### [`robin utils`](utils.md)
Models, ClinVar, `sequencing-files`, `mgmt`, and more.
</div>

<div class="robin-feature-card" markdown>
### [`robin password`](password.md)
Set or replace the NiceGUI login password.
</div>

</div>

---

## Disclaimer (`I agree`)

Commands that start processing or change protected state usually show a **research-use disclaimer** and require typing **`I agree`** exactly.

---

## Environment variables (selected)

| Variable | Effect |
|----------|--------|
| `GITHUB_TOKEN` | Private GitHub assets for `robin utils update-models`. |
| `ROBIN_PROCESS_LARGE_BAMS` | Warns that large-BAM mode must not run alongside live sequencing. |
| `LJ_BAM_THREADS` | Optional BAM threading ([README — Performance](https://github.com/LooseLab/ROBIN/blob/main/README.md#performance)). |

---

## Startup behaviour

**`robin workflow`** runs model checks, optional reference validation, **`I agree`**, and (if the GUI is on) **GUI password** prompts — see **[What happens at startup](../getting-started/startup.md)**.

---

## Web interface

Browser guide: **[Using ROBIN](../using-robin/index.md)**.

---

## See also

- [README — Command reference](https://github.com/LooseLab/ROBIN/blob/main/README.md#command-reference)
