# MinKNOW configuration

How to configure MinKNOW for real-time BAM output compatible with ROBIN.

## Overview

ROBIN expects:

- Real-time **basecalling** (with **methylation** output as required by your ROBIN build)
- **Alignment** to the same reference you pass to ROBIN
- **BAM** output with sensible rollover (read-count–based, not hourly-only)
- Clear **sample** (and optionally **experiment**) IDs for output organisation

## Basecalling

Use **high-accuracy** basecalling with **5hmC / 5mc** methylation calling where ROBIN requires it. **Do not** disable methylation for off-target reads if your pipeline uses those signals.

Example basecalling config (verify against your MinKNOW / kit release):

- `dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg`

## Adaptive sampling

If you use adaptive sampling:

- Use the **BED** file specified for your ROBIN / panel build (e.g. `panel_11092024_5kb_pad.bed` for classic CNS panels — confirm path in your installed resources).
- Use the **same reference** as ROBIN and alignment.
- Mode typically **enrich** (per your assay design).

## Alignment

Produce **aligned BAMs** against the **same** reference genome files used in ROBIN (`--reference`).

## BAM output

- Prefer **read-count** rollover, not **time-only** hourly chunks (unsupported for real-time ROBIN processing).
- Keep each BAM **under the supported read count** (see [README](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings)); a common target is **~20,000–50,000 reads per file** depending on workflow.
- You may disable **POD5** and **FASTQ** if you only need BAM.

## Sample and experiment IDs

- Set a **unique sample ID** per library — ROBIN uses this for output folders and to join analysis to the correct run.
- If you use ROBIN’s **[Sample ID generator](../using-robin/pages-and-routes.md#sample-id-generator)** (MD5 from test ID and optional patient fields), use that **generated sample ID** as the **MinKNOW sample ID** when you start the run. ROBIN will then **link** MinKNOW outputs, the on-disk manifest, and the web app **without** manual renaming.
- **Experiment ID** can group runs; you can point ROBIN at a subset of output directories if helpful.
- Set **run duration** appropriate to your clinical or research protocol (e.g. 24 h for many CNS assays).

## Run workflow

1. Apply the settings above in MinKNOW.  
2. Start the run.  
3. Start ROBIN with the directory MinKNOW writes BAMs into (`robin workflow …` for Little John, or legacy `robin …` for Original ROBIN).  
4. Monitor progress in the CLI / NiceGUI (Little John) or legacy GUI (Original ROBIN).

## Troubleshooting

- **API errors** (`minknow_api`): install the Python `minknow-api` version that matches your MinKNOW install — see [Installation → MinKNOW API](installation.md).  
- **Huge BAMs / missed processing**: reduce reads per file and avoid time-only rollover.  
- **Reference mismatch**: alignment reference must match ROBIN’s `--reference`.

## Next steps

- [Quickstart](quickstart.md)  
- [MinKNOW documentation](https://nanoporetech.com/minknow)  
- [Adaptive sampling](https://nanoporetech.com/adaptive-sampling)
