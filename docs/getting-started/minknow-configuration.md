# MinKNOW configuration

!!! abstract "What this page covers"
    Configure **MinKNOW** so it streams **aligned BAMs** that ROBIN can consume: basecalling (with **methylation** if required), **same reference** as ROBIN, **read-count–based** BAM rollover, and clear **sample IDs**.  
    Then start ROBIN on the folder MinKNOW writes BAMs into — [Quickstart](quickstart.md).

---

## Optional: programmatic MinKNOW integration

ROBIN can connect to MinKNOW over the API to show live sequencer status and auto-watch BAM output directories. This requires the optional Python client:

```bash
pip install 'robin[minknow]'
```

The **`minknow_api`** package **minor version must match MinKNOW Core** on your instrument (check **Host Settings → About**). Examples:

| MinKNOW Core | Install |
|--------------|---------|
| 6.8.x | `pip install 'minknow_api>=6.8.0,<6.9.0'` |
| 6.10.x | `pip install 'minknow_api>=6.10.0,<6.11.0'` |

CLI: `robin minknow status --host <sequencer>`. GUI: **Sequencer (MinKNOW)** card on the workflow and live-data pages. See [MinKNOW integration plan](../development/minknow-integration.md) for auth, auto-watch, and env vars (`MINKNOW_HOST`, `MINKNOW_AUTO_WATCH`, etc.).

---

## What ROBIN needs

| Area | Requirement |
|------|-------------|
| **Basecalling** | Real-time; enable **5mC / 5hmC modified-base calling in CpG contexts only** where your ROBIN build needs methylation. Do **not** use all-context calling. |
| **Alignment** | Same **reference** you pass to **`robin workflow --reference`**. |
| **BAM output** | **Read-count** rollover — **not** time-only hourly chunks. Keep each file under the [supported read count](https://github.com/LooseLab/ROBIN/blob/main/README.md#bam-read-limit-and-minknow-settings); **~20k–50k reads per file** is a common range. |
| **Sample ID** | Unique per library; match [generated IDs](../using-robin/pages-and-routes.md#sample-id-generator) if you use the Sample ID helper. |

---

## Basecalling

Use **high-accuracy (HAC)** basecalling with **5mC / 5hmC modified-base calling in CpG contexts only** where ROBIN requires methylation. Do **not** select an all-context 5mC / 5hmC model. Do not disable CpG methylation calling for off-target reads if your pipeline uses those signals.

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
| **`minknow_api` errors** | Install **`robin[minknow]`** and a **`minknow_api`** version whose **minor version matches MinKNOW Core** (see [Optional: programmatic MinKNOW integration](#optional-programmatic-minknow-integration) above). |
| **Huge BAMs / missed files** | Reduce reads per file; avoid time-only rollover. |
| **Reference mismatch** | Alignment reference must match **`--reference`**. |

---

## Next steps

- [Quickstart](quickstart.md)  
- [MinKNOW product docs](https://nanoporetech.com/minknow)  
- [Adaptive sampling](https://nanoporetech.com/adaptive-sampling)  
