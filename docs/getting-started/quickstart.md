# ROBIN Quickstart

!!! abstract "What this page covers"
    Run **`robin workflow`** after [installation](installation.md): what ROBIN expects from **MinKNOW/BAMs**, how to stage **reference + panel BEDs**, example commands, and **panels**.  
    For **disclaimer, GUI password, and startup order**, see **[What happens at startup](startup.md)**. For the **web UI**, see **[Using ROBIN](../using-robin/index.md)**. For every flag, see **[CLI reference](../cli/index.md)** and the **[README](https://github.com/LooseLab/ROBIN/blob/main/README.md)**.

---

## What ROBIN expects from sequencing

ROBIN consumes **aligned BAMs** from Oxford Nanopore (usually written in real time by MinKNOW):

| Expectation | Notes |
|-------------|--------|
| **Basecalling** | **HAC** is enough; SUP not required. |
| **Methylation** | Enable **5hmC / 5mC** in MinKNOW if your analyses need it. |
| **Alignment** | Done **in MinKNOW** — ROBIN does not realign. |
| **BAM chunks** | **Read-count–based** rollover. Each BAM **≤ 50,000 reads**; **~50k per file** is a good target. **Do not** use **time-only** (e.g. hourly) rollover — see [README — BAM read limit](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings). |
| **POD5 / FASTQ** | Optional; you can disable if BAM alone is enough. |

!!! tip "Memory on smaller machines"
    On **≤ 64 GB RAM**, restart between long runs or after moving the flow cell. Dorado can retain GPU/host memory; restarting Dorado or the instrument after a run reduces OOM risk.

---

## Reference and panel BEDs (one command) {#reference-and-panel-beds-one-command}

Use **`robin utils sequencing-files`** to gather **panel BED** + **GRCh38 reference** in one folder for MinKNOW and for **`robin workflow --reference`**. Run after [installation](installation.md) when you need a single consistent reference path.

```bash
robin utils sequencing-files --panel rCNS2 --output-dir ~/references/robin_ref
```

| Option | Purpose |
|--------|---------|
| `-p` / `--panel` | **Required.** Same names as `--target-panel` (`rCNS2`, `AML`, `PanCan`, … — run `robin utils sequencing-files --help` on your install). |
| `-r` / `--reference` | **Reference FASTA:** HTTPS URL or local `.fa` / `.fa.gz`. If omitted, ROBIN may download the default **GRCh38 no-alt** set (large). |
| `-o` / `--output-dir` | Output folder (default: `./reference_files`). |
| `-y` / `--yes` | Skip confirmation (for scripts). |

Use the **same** reference file for **MinKNOW alignment** and **`robin workflow --reference`**.

### Other `robin utils` commands

| Command | Purpose |
|---------|---------|
| `robin utils update-models` | Models / classification assets ([Installation](installation.md)). |
| `robin utils update-clinvar` | ClinVar resources for annotation. |
| `robin utils mgmt` | Summarise **MGMT** CpG methylation from existing `mgmt_sorted.bam` outputs. |

Run **`robin utils --help`** for the full list.

---

## Run a workflow

Typical invocation:

```bash
robin workflow <data_folder> --work-dir <output_folder> \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center <center_id>
```

| Argument | Meaning |
|----------|---------|
| `<data_folder>` | Directory watched for incoming BAMs |
| `--work-dir` | Root for outputs |
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

Point `--reference` at the **same GRCh38** file you use in MinKNOW—for example from your [`sequencing-files`](#reference-and-panel-beds-one-command) output.

---

## Commands you’ll use often

### `list-job-types`

```bash
robin list-job-types
```

Examples include preprocessing, bed_conversion, mgmt, cnv, target, fusion, sturgeon, nanodx, pannanodx, random_forest (see live output for your version).

### `workflow`

```bash
robin workflow /path/to/directory -w "<workflow_plan>" [OPTIONS]
```

**Required:** `-w` / `--workflow`, and `--center`.

**Useful:** `--work-dir`, `--reference`, `--verbose`, `--log-level`, `--no-process-existing`, `--deduplicate-jobs`, `--use-ray` / `--no-use-ray`, `--with-gui` / `--no-gui` — full list: **`robin workflow --help`**.

---

## Panel management

Built-in panels include **rCNS2**, **AML**, **PanCan**. Custom panels are registered from BED (≥ four columns: chr, start, end, gene name(s)).

```bash
robin list-panels
robin add-panel /path/to/panel.bed MyCustomPanel
robin add-panel /path/to/panel.bed MyCustomPanel --validate-only
robin remove-panel MyCustomPanel
robin remove-panel MyCustomPanel --force
```

Reserved names: **`rCNS2`**, **`AML`**, **`PanCan`**.

---

## Behaviour and limits

- **CNV** — Heuristic; **review visually** before clinical use.  
- **Stop with Ctrl+C** — Graceful shutdown is attempted but not guaranteed.  
- **Issues** — [GitHub issues](https://github.com/LooseLab/ROBIN/issues).  

### Performance (brief)

Batched processing on heavy paths; optional **`LJ_BAM_THREADS`**; non-blocking GUI when NiceGUI is enabled.

### License

Research use; see **LICENSE** in the repository. ROBIN integrates tools such as Sturgeon, Rapid-CNS2, Readfish, and others—see the repo for attribution.

---

## Next steps

<div class="grid" markdown>

<div class="grid-item" markdown>
### Wet lab
[Library preparation](library-preparation.md)
</div>

<div class="grid-item" markdown>
### Instrument
[MinKNOW configuration](minknow-configuration.md)
</div>

<div class="grid-item" markdown>
### Deep dive
[README — Usage](https://github.com/LooseLab/ROBIN/blob/main/README.md#usage)
</div>

</div>
