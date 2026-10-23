# Experimental live SV calling and replay

ROBIN can accumulate structural-variant evidence from aligned BAM batches using
the optional `sv` workflow stage. The sample's Live Data page provides a read-only
SV viewer, mapping metrics, an experimental evidence score and linked CNV highlights.
This is a research prototype: evidence scores are not probabilities, somatic
classification or validated diagnostic calls.

## Enable the native worker

Add `sv` to your existing workflow selection, for example
`cnv,sturgeon,nanodx,mgmt,fusion,sv`. The native Ray coordinator schedules SV work
after preprocessing in its shared analysis pool. Existing workflows do not gain
SV analysis unless requested.

Results live under `<work-dir>/<sample>/sv_live/`:

- `state.sqlite`: persistent segments, oriented junction clusters, unique-read
  evidence, masks, local-anchor quality and mapping summaries.
- `updates.jsonl`: changed calls, metrics and versioned evidence scores.
- `batches.jsonl`: input paths and service timings.

Writers use a per-sample file lock and SQLite transactions. Reprocessing the same
read does not inflate support. Primary/supplementary records can arrive in separate
batches; changed reads are reconsidered. Robust publication-feed recovery after
arbitrary interruption has not been established. Keep one reference dictionary
and one biological sample per checkpoint.

## Replay archived batches through the full pipeline

Prepare a manifest from the original sequencing summary and BAM batch directories:

```bash
python -m robin.sv_replay.benchmark prepare \
  --summary /path/to/sequencing_summary.txt \
  --bam-dir /path/to/bam_pass --bam-dir /path/to/bam_fail \
  --duration 86400 --output /path/to/manifest.csv

python -m robin.sv_replay.full_run \
  --manifest /path/to/manifest.csv \
  --reference /path/to/matching_reference.fa \
  --output /path/to/new_replay --speed 1
```

The reference FASTA is required explicitly. Manifest columns are `bam,available_at`.
Availability is estimated from each original BAM's latest read completion; it is
not the historical basecalling/alignment/publication time. Initial idle time is
preserved. `--speed 1` uses original pacing; faster replay increases contention.

The replay uses ROBIN's native watcher, preprocessing, CNV, bed conversion,
Sturgeon, NanoDX, MGMT, fusion and SV stages. Required classifier assets must be
installed separately through ROBIN's existing asset workflow. The output must be
new and empty. Source BAMs are published as symlinks and are not modified.

Open a separate GUI on the replay's work directory:

```bash
python -m robin.gui_launcher /path/to/new_replay/robin \
  --host 127.0.0.1 --port 8091 --no-show
```

Use the existing ROBIN account and approval controls. The sample directory appears
only after the first batch is processed. Restarting the GUI does not restart the
analysis. Completed runs include a summary and final TSV exports. Full-run replay
resume is not implemented.

## Standalone offline evidence replay

```bash
python -m robin.sv_replay \
  --bam /path/to/aligned.bam \
  --read-times /path/to/read_times.csv \
  --output /path/to/new_evidence_replay
```

After installation, `robin-sv-replay` is equivalent. The timing CSV contains
`read_id,available_at`; every retained read must have a finite non-negative time.
The offline replay stages all input before consuming reads in time order. It runs
as fast as possible, not at wall-clock pacing. Without timing, repeated `--bam`
arguments use synthetic batch timing; a coordinate-sorted full BAM is not a live
arrival-order replay.

Outputs are `replay.sqlite`, `calls.tsv`, `evidence.tsv`, `history.tsv` and `run.json`.
Reissuing an unchanged offline command resumes committed work. Input/configuration
identity is bound to the checkpoint; a changed policy requires a fresh directory.
Mapping/error summaries and score v2 are supplied by the incremental worker;
standalone offline output does not automatically provide those live metrics.

## Evidence and annotations

Discovery requires MAPQ >=20, a >=100-base segment anchor, and SV size >=50 bp.
CIGAR I/D and compatible query-adjacent split alignments generate signatures;
query gap/overlap is limited to 100 bp. Junctions cluster with matching oriented
ends within 100 bp of a fixed seed. Equivalent DEL and BND representations are
not merged automatically. Parent/duplex relationships are not reconciled.

Breakpoints use zero-based interbase coordinates. Retained-side orientations
classify same-chromosome junctions as `deletion_like` (+/-), `duplication_like`
(-/+), or `inversion_like` (equal signs). These are junction topologies, not fully
reconstructed events. Inversions, complex amplifications and copy counts require
additional evidence.

Default GRCh38 annotations are:

- `resources/cenSatRegions.bed`: the existing broad centromere/satellite context.
- `resources/sv_grch38_tandem_repeats.bed`: tandem repeats from Severus; source,
  checksum and upstream license accompany the file.
- A 1 Mb chromosome-end margin; a special mitochondrial-reference-origin rule.

Only the retained breakpoint flank is annotated, not the entire event span.
Tandem-repeat and centromere annotations are context, not evidence of false calls.
They do not reduce score v2. Status remains separate: fewer than 3 IDs is
`candidate`; >=3 is `provisional` when unmasked or `masked_candidate` with masks.
Mitochondrial origin junctions use `reference_boundary`.

The standalone CLI accepts `--tandem-repeats`, `--centromeres`, extra `--mask-bed`
annotations and configurable end/origin margins. Native defaults are GRCh38;
do not use those masks for another assembly without adapting the worker.

## Mapping metrics and evidence score v2

The worker records MAPQ on each junction side and a whole-segment substitution
estimate `(NM - CIGAR I/D bases) / aligned M/=/X bases`. Missing or impossible NM
is unknown. Actual alignment records are preferred to approximate SA descriptions
at matching chromosome/start/strand. The estimate is not breakpoint-local identity.

Summaries include measured ID/endpoint counts, MAPQ >=30 on both sides, mean
minimum MAPQ, mean substitution rate, and population SD of both breakpoints and
event size. Interchromosomal event size is unknown. SD=0 for one observation does
not establish consistency. The strict local-anchor metric counts consecutive M/=
bases beside each end; both must be >=100 bp. Small nanopore errors can interrupt
these anchors, and M does not guarantee exact identity.

The score is implemented in `sv_replay/scoring.py`:

| Component | Maximum | Rule |
|---|---:|---|
| Support | 30 | Logarithmic increase, saturating at 10 IDs |
| Mapping | 25 | Fraction of total IDs with both MAPQs >=30 |
| Local anchors | 10 | Fraction meeting the strict anchor requirement |
| Substitutions | 20 | Decreases to zero at 8%; weighted by measurement coverage |
| Consistency | 15 | Requires >=3 measured IDs; depends on position/size dispersion |

One/two IDs cap at 25/39. Incomplete mapping/error measurements cap at 69, as do
<3 strong-mapping IDs, <50% strong mapping, substitution burden >5%, coordinate
SD >50 bp, excessive size dispersion, or other custom/end-region flags.
Reference-boundary calls cap at 5. Centromere and tandem-repeat tags have no score
penalty. Low/medium/high export categories correspond to <40 / 40–69 / >=70.
These thresholds are heuristic and need independent truth-set validation.

Native worker updates collect these metrics directly. For legacy checkpoints,
`python -m robin.sv_replay.backfill_metrics /path/to/sample/sv_live --follow` can
recover them from archived BAM paths recorded in `batches.jsonl`. This auxiliary
pass rereads BAMs, writes derived tables only, and follows arrivals every 30–60 s.
It does not alter calls, masks, statuses or supporting IDs.

## Viewer controls

The viewer defaults to >=3 supporting IDs and hides measured mean minimum MAPQ <55.
The checkbox **Show low mapping quality calls (MAPQ<55)** includes them. Unknown
MAPQ stays visible. Sorting is numeric and applies before pagination. The green
score bar and information popup expose the scoring factors; detailed metrics are
available in the popup/download rather than expanding every column.

Clicking a call switches the CNV chromosome, scrolls to the profile and highlights
its reference interval in translucent purple. Spans <10 kb use a midpoint line;
interchromosomal calls get separate endpoint bands. Highlights hold 1.5 seconds
then fade over 1.5 seconds, using a DOM overlay rather than chart redraws.
**Supporting reads** opens up to 100 IDs with per-read measurements. Page downloads
use ROBIN's existing export authorization. CNV data and analysis are not modified
by selecting an SV.

## Verification and limitations

```bash
python -m unittest discover -s tests -p 'test_sv*.py'
python -m unittest discover -s tests -p 'test_cnv_sv_refresh.py'
```

Use an installed ROBIN environment with its optional GUI/native dependencies for
integration tests. Synthetic tests cover extraction, retries, masks, metrics,
scoring, filtering, CNV highlight coordinates and transient-file handling. They
are not a biological validation or full-pipeline performance guarantee. Repeat
normalization, alternative-placement searches, somatic classification and general
complex-event reconstruction are not implemented. Keep sequencing data, sessions,
account databases, replay outputs and downloaded model assets out of Git.
