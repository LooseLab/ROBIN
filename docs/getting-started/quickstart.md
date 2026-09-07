# ROBIN Quickstart

Use this page to get from an installed ROBIN environment to a running analysis with the web interface open.

!!! note
    If ROBIN is not installed yet, start with [Installation](installation.md).

## 1. Prepare the reference and panel

Create a folder containing the reference FASTA and the BED files required for your sequencing setup:

```bash
robin utils sequencing-files \
  --panel rCNS2 \
  --output-dir ~/references/robin_ref
```

By default this stages the built-in `rCNS2` panel and downloads the supported GRCh38 no-alt reference. If you already have a compatible GRCh38 FASTA, provide it with `--reference`.

Use the **same reference FASTA** for MinKNOW alignment and for `robin workflow --reference`.

For all options, see [`robin utils sequencing-files`](../cli/utils.md#robin-utils-sequencing-files).

## 2. Configure MinKNOW

ROBIN expects aligned Oxford Nanopore BAM files produced during sequencing.

For the standard real-time workflow:

- use HAC basecalling or better;
- enable **5mC/5hmC calling in CpG context** when methylation analyses are required;
- align reads in MinKNOW using the same reference supplied to ROBIN;
- configure BAM output by **read count**, not time;
- keep each BAM at **50,000 reads or fewer**.

!!! warning "BAM rollover matters"
    Time-based BAM rollover can produce files that are too large for ROBIN's supported real-time workflow. Configure MinKNOW to rotate BAMs by read count; approximately 50,000 reads per BAM is recommended.

For the complete instrument setup, see [MinKNOW configuration](minknow-configuration.md).

## 3. Start ROBIN

Assuming MinKNOW is writing BAMs to `~/data/bam_files` and the reference is at `~/references/robin_ref/hg38.fa`:

```bash
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/robin_ref/hg38.fa \
  --center Sherwood \
  --target-panel rCNS2
```

Replace the paths and `--center` value for your site.

ROBIN automatically adds preprocessing and other required upstream steps where appropriate. To see the analysis names available in your installed version, run:

```bash
robin list-job-types
```

For detailed workflow options, see [`robin workflow`](../cli/workflow.md) and [Job types](../cli/jobs.md).

## 4. Complete startup

On first use, ROBIN may ask you to:

1. accept the research-use disclaimer;
2. create the initial GUI administrator password.

With the GUI enabled, ROBIN prints the web address it is serving. Open that address in a browser and sign in.

For the exact startup sequence, authentication behaviour, Ray/threading mode and shutdown behaviour, see [Starting ROBIN](startup.md).

## 5. Confirm that the run is working

Once BAMs appear in the watched directory, confirm that:

- ROBIN detects the sample;
- preprocessing jobs complete;
- the sample appears in the web interface;
- enabled analyses begin producing results as data accumulate.

Then continue with:

- [Using ROBIN](../using-robin/index.md) — navigation, samples and the web interface;
- [Reading your results](../using-robin/sample-results.md) — interpreting the GUI;
- [Analysis pipelines](../analyses/index.md) — what each analysis does;
- [Troubleshooting](../using-robin/troubleshooting.md) — common operational problems.

## Common variations

### Run only selected analyses

You do not need to enable every analysis. For example:

```bash
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/robin_ref/hg38.fa \
  --center Sherwood \
  --target-panel rCNS2
```

### Use a different panel

Built-in and custom panel management is documented under [Panel commands](../cli/panels.md).

### Run without the browser interface

```bash
robin workflow ... --no-gui
```

See the [workflow command reference](../cli/workflow.md) for other execution options.
