# `robin list-job-types` and job model

`robin list-job-types` prints the **job types** ROBIN knows about, grouped by **queue**. Understanding this helps you build `-w` / `--workflow` strings for [`robin workflow`](workflow.md).

## Command

```bash
robin list-job-types
```

You must complete the **disclaimer** (`I agree`) unless consent has already been recorded for an active administrator.

## Queues and job types

The `release_candidate2` CLI registers the following workflow handlers:

| Queue | Job types | Role |
| --- | --- | --- |
| **Preprocessing** | `preprocessing` | Validate each BAM and extract sample/run metadata. |
| **BED conversion** | `bed_conversion` | Extract modified-base data and maintain classifier input. |
| **MGMT** | `mgmt` | MGMT promoter methylation analysis. |
| **CNV** | `cnv` | Genome-wide copy-number analysis. |
| **Target** | `target` | Target coverage and downstream variant-analysis preparation. |
| **Fusion / structural** | `fusion`, `itd` | Fusion/rearrangement and ITD/insertion analysis. |
| **Classification** | `sturgeon`, `nanodx`, `pannanodx` | Methylation classifiers. |
| **Slow** | `random_forest`, `marlin`, `lamprey`, `tucan` | Heavier/specialised classifiers and research analyses. |

The older queue-qualified workflow syntax may use legacy queue names internally; the simplified job list above is preferable for normal use.

## Automatic steps

When you use the simplified workflow format (`-w mgmt,sturgeon`, for example), ROBIN constructs the required pipeline around the requested analyses.

In particular, classifier jobs that consume the methylation parquet require `bed_conversion` upstream. Preprocessing is the entry point for new BAM files and establishes metadata used by downstream jobs.

You normally do not need to manually construct the legacy queue-prefixed form.

## Valid job type names

The CLI handler configuration in `release_candidate2` includes:

`preprocessing`, `bed_conversion`, `mgmt`, `cnv`, `target`, `fusion`, `itd`, `sturgeon`, `nanodx`, `pannanodx`, `random_forest`, `marlin`, `lamprey`, `tucan`

Use `robin list-job-types` as the authoritative runtime list for your installed checkout.

## Examples

```bash
# Minimal classifier run
robin workflow /data/bams \
  -w sturgeon \
  --center Demo \
  --target-panel rCNS2 \
  -d /out \
  --reference /ref/hg38.fa

# Broad CNS analysis
robin workflow /data/bams \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --center Sherwood \
  --target-panel rCNS2 \
  -d ~/results \
  --reference ~/references/hg38.fa

# Add ITD analysis where the active panel/configuration defines applicable hotspots
robin workflow /data/bams \
  -w target,itd \
  --center Demo \
  --target-panel AML \
  -d /out \
  --reference /ref/hg38.fa
```

## ITD notes

`itd` shares the structural/fusion scheduling path but has its own handler. It requires an active target panel and resolves scan windows from the ITD hotspot configuration. If no applicable hotspots overlap the active panel, ROBIN records an empty result rather than scanning arbitrary genomic regions.

See [Structural events: fusions and ITDs](../analyses/structural-events.md).

## MARLIN notes

- Install the optional extra: `pip install 'robin[marlin]'` where required by your checkout.
- Model/runtime requirements are separate from the core ROBIN workflow.
- Treat MARLIN output according to the limitations and licensing of the upstream model.

## Lamprey notes

Lamprey support is intended for research/evaluation and has upstream licensing/model requirements. Install and use it only where those requirements are satisfied.

## Experimental/specialised jobs

`marlin`, `lamprey`, and `tucan` are specialised paths and may have additional model, dependency, licensing, or validation requirements. Their presence in `list-job-types` does not imply that all required external resources are installed.

## Related

- [`robin workflow`](workflow.md)
- [Quickstart](../getting-started/quickstart.md)
- [Analysis pipelines](../analyses/index.md)
- [ROBIN architecture](../architecture/overview.md)
