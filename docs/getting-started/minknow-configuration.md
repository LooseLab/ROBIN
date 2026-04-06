# MinKNOW configuration

!!! abstract "What this page covers"
    Configure **MinKNOW** so it streams **aligned BAMs** that ROBIN can consume: basecalling (with **methylation** if required), **same reference** as ROBIN, **read-count–based** BAM rollover, and clear **sample IDs**.  
    Then start ROBIN on the folder MinKNOW writes BAMs into — [Quickstart](quickstart.md).

---

## What ROBIN needs

| Area | Requirement |
|------|-------------|
| **Basecalling** | Real-time; **5hmC / 5mC** where your ROBIN build needs methylation. |
| **Alignment** | Same **reference** you pass to **`robin workflow --reference`**. |
| **BAM output** | **Read-count** rollover — **not** time-only hourly chunks. Keep each file under the [supported read count](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings); **~20k–50k reads per file** is a common range. |
| **Sample ID** | Unique per library; match [generated IDs](../using-robin/pages-and-routes.md#sample-id-generator) if you use the Sample ID helper. |

---

## Basecalling

Use **high-accuracy** basecalling with **5hmC / 5mc** methylation calling where ROBIN requires it. Do **not** turn off methylation for off-target reads if your pipeline uses those signals.

Example config name (verify against your MinKNOW / kit release):

- `dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg`

---

## Adaptive sampling

If you use adaptive sampling:

- Use the **BED** for your ROBIN / panel build (e.g. panel BED from resources — confirm path for your install).  
- Use the **same reference** as ROBIN and MinKNOW alignment.  
- Mode is typically **enrich** per your assay design.  

---

## Alignment

Produce **aligned BAMs** against the **same** FASTA files as **`robin workflow --reference`**.

---

## BAM output

- Prefer **read-count** rollover, not **time-only** chunks.  
- Keep each BAM within the **supported read limit** (see [README](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings)).  
- You may disable **POD5** and **FASTQ** if only BAM is needed.  

---

## Sample and experiment IDs

- **Sample ID** — unique per library; ROBIN uses it for output folders and tracking.  
- **Sample ID generator** — if you use ROBIN’s [MD5-based ID](../using-robin/pages-and-routes.md#sample-id-generator), enter that value as the **MinKNOW sample ID** so outputs line up automatically.  
- **Experiment ID** — optional grouping of runs.  
- **Run duration** — set per your clinical or research protocol (e.g. 24 h for many CNS assays).  

---

## Suggested workflow

1. Apply the settings above in MinKNOW.  
2. Start the run.  
3. Start **`robin workflow …`** pointing at the directory MinKNOW writes BAMs into.  
4. Monitor in the CLI / NiceGUI ([Using ROBIN](../using-robin/index.md)).  

---

## Troubleshooting

| Symptom | What to check |
|---------|----------------|
| **`minknow_api` errors** | Install a **`minknow-api`** version that matches your MinKNOW install ([Installation](installation.md) if documented). |
| **Huge BAMs / missed files** | Reduce reads per file; avoid time-only rollover. |
| **Reference mismatch** | Alignment reference must match **`--reference`**. |

---

## Next steps

- [Quickstart](quickstart.md)  
- [MinKNOW product docs](https://nanoporetech.com/minknow)  
- [Adaptive sampling](https://nanoporetech.com/adaptive-sampling)  
