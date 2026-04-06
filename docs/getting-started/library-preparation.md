# Library preparation

!!! abstract "What this page covers"
    **Oxford Nanopore** library choices for work later analysed with ROBIN—read length targets and links to the **ROBIN protocol** vs standard **LSK114**.  
    Next: [MinKNOW configuration](minknow-configuration.md), then [Quickstart](quickstart.md).

---

## Choose a protocol

ROBIN is designed for **Oxford Nanopore** libraries. For CNS tumour work you can use either:

| Option | Where to read more |
|--------|---------------------|
| **ROBIN-specific protocol** | [ROBIN: unified nanopore assay (protocols.io)](https://www.protocols.io/view/robin-a-unified-nanopore-based-sequencing-assay-in-bp2l6xepklqe/v3) — fragmentation, ligation, QC, troubleshooting. |
| **Standard ONT LSK114** | [LSK114 product page](https://store.nanoporetech.com/uk/ligation-sequencing-kit-v14.html) and [library preparation hub](https://nanoporetech.com/library-preparation). |

---

## Read length

**Adaptive sampling** works best with read lengths around **8–12 kb**—follow your chosen protocol to stay in that range.

---

## DNA input

- High-quality DNA (e.g. A260/A280 ~1.8–2.0)  
- Minimum ~**100 ng**; often **500 ng–1 µg** is optimal  
- Elute in TE or nuclease-free water  

---

## Next steps

1. [MinKNOW configuration](minknow-configuration.md) — basecalling, alignment, BAM rollover.  
2. [Quickstart](quickstart.md) — run ROBIN on the output directory.  

---

## Resources

- [ROBIN protocol (protocols.io)](https://www.protocols.io/view/robin-a-unified-nanopore-based-sequencing-assay-in-bp2l6xepklqe/v3)  
- [ONT library preparation](https://nanoporetech.com/library-preparation)  
- [SQK-RAD114](https://store.nanoporetech.com/sqk-rad114.html) (where applicable)  
