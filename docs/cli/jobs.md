# `robin list-job-types` and job model

`robin list-job-types` prints the **job types** ROBIN knows about, grouped by **queue**. Understanding this helps you build `-w` / `--workflow` strings for [`robin workflow`](workflow.md).

## Command

```bash
robin list-job-types
```

You must complete the **disclaimer** (`I agree`).

## Queues and job types

The orchestration layer assigns each job type to a queue (simplified names here; internal Ray queue names may differ slightly):

| Queue (concept) | Job types | Role |
|-----------------|-----------|------|
| **Preprocessing** | `preprocessing` | Read BAM headers/metadata; entry point for each new file. |
| **BED conversion** | `bed_conversion` | Prepare inputs for classifiers that need BED-level views. |
| **Analysis** | `mgmt`, `cnv`, `target`, `fusion` | Methylation (MGMT), copy number, targeted variant/fusion panels. |
| **Classification** | `sturgeon`, `nanodx`, `pannanodx` | Methylation / expression classifiers. |
| **Slow** | `random_forest`, `marlin`, `lamprey` | Heavier models (RF; MARLIN TF; Lamprey ONNX research-only). |

## Automatic steps

When you use the **simplified** workflow format (`-w mgmt,sturgeon`, …):

1. **`preprocessing`** is prepended if you did not list it.
2. **`bed_conversion`** is inserted when any of **`sturgeon`**, **`nanodx`**, **`pannanodx`**, **`random_forest`**, **`marlin`**, or **`lamprey`** appear — those jobs expect BED conversion upstream.

You do not need to list `bed_conversion` manually for those classifiers unless you are hand-editing **legacy** queue-prefixed pipelines.

## Valid job type names

The CLI accepts only these **job** identifiers in workflow strings:

`preprocessing`, `bed_conversion`, `mgmt`, `cnv`, `target`, `fusion`, `sturgeon`, `nanodx`, `pannanodx`, `random_forest`, `marlin`, `lamprey`

Unknown names produce warnings and are skipped.

## Examples

```bash
# Minimal classifier run (preprocessing + bed_conversion added as needed)
robin workflow /data/bams -w sturgeon --center Demo --target-panel rCNS2 -d /out --reference /ref/hg38.fa

# Full stack (typical)
robin workflow /data/bams \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest,marlin,lamprey \
  --center Sherwood \
  --target-panel rCNS2 \
  -d ~/results \
  --reference ~/references/hg38.fa
```

## MARLIN notes

- Install the optional extra: `pip install 'robin[marlin]'` (pulls TensorFlow ≥2.16 and `tf-keras` for Keras-2 HDF5 loading; required for Python 3.12+).
- On first MARLIN job, ROBIN downloads `marlin_v1.model.hdf5` (~1.1 GiB) from Zenodo into `~/.cache/robin/marlin/` (override with `ROBIN_MARLIN_MODEL_PATH` / `ROBIN_MARLIN_CACHE_DIR`).
- Default probe genome build is **hg38** (also supports `hg19` / `t2t` via job metadata `marlin_genome_build`).

## Lamprey notes (research / evaluation only)

**Lamprey is not for clinical care, diagnosis, or medical decision-making.** Its upstream license restricts use to internal non-commercial research and evaluation (Oncode / Cyclomics / UMCU). Contact `software@cyclomics.com` for clinical/commercial licensing.

- Install Lamprey **separately** (ROBIN does not vendor it), e.g. `pip install git+ssh://git@github.com/princessmaximacenter/lamprey.git`, plus `pip install 'robin[lamprey]'`.
- Acknowledge the research terms before first model download: `export ROBIN_LAMPREY_RESEARCH_ACK=1`.
- On first run, ROBIN downloads the HuggingFace model (`tachterberg/Lamprey`, ~7 GiB) into `~/.cache/robin/lamprey/` (override with `ROBIN_LAMPREY_MODEL_PATH` / `ROBIN_LAMPREY_CACHE_DIR`).
- ROBIN Lamprey support is **hg38 only**. Confidence tiers match Sturgeon (high ≥95%, medium ≥80%).

## Related

- [`robin workflow`](workflow.md)  
- [Quickstart](../getting-started/quickstart.md)  
