# ROBIN

![ROBIN Logo](images/ROBIN_logo_small.png)

**Rapid nanopOre Brain intraoperatIve classificatioN**

!!! warning "Research use only"
    ROBIN is for **research use** at this time. The technology is under active development and validation.

!!! abstract "What ROBIN is"
    ROBIN is a real-time analysis and web-monitoring platform for Oxford Nanopore sequencing of CNS tumours. It combines live BAM processing with methylation classification, copy-number analysis, MGMT analysis, targeted analysis and structural-event detection.

## Get started

<div class="grid" markdown>

<div class="robin-feature-card" markdown>
### Install ROBIN
[Installation](getting-started/installation.md) covers the conda environment, package installation and required model assets.
</div>

<div class="robin-feature-card" markdown>
### Run your first workflow
[Quickstart](getting-started/quickstart.md) takes you from an installed environment to a running analysis and web interface.
</div>

<div class="robin-feature-card" markdown>
### Configure sequencing
[MinKNOW configuration](getting-started/minknow-configuration.md) covers alignment, modified-base calling and BAM rollover settings.
</div>

</div>

## Explore ROBIN

<div class="grid" markdown>

<div class="robin-feature-card" markdown>
### Use the web interface
[Using ROBIN](using-robin/index.md) covers sign-in, navigation, samples, results and troubleshooting.
</div>

<div class="robin-feature-card" markdown>
### Understand the analyses
[Analysis pipelines](analyses/index.md) explains the purpose, inputs and interpretation of ROBIN's analysis workflows.
</div>

<div class="robin-feature-card" markdown>
### Command-line reference
[CLI reference](cli/index.md) documents `robin workflow`, job types, panels and utilities.
</div>

</div>

---

## Why rapid molecular analysis matters

Brain and CNS tumours encompass many biologically distinct entities, and molecular information increasingly contributes to classification and research workflows. Conventional diagnostic pathways may take days to weeks to assemble the complete molecular picture.

ROBIN was developed to explore whether nanopore sequencing and real-time analysis can shorten that interval substantially. The approach is described in our [Neuro-Oncology publication](https://academic.oup.com/neuro-oncology/advance-article/doi/10.1093/neuonc/noaf103/8139084?searchresult=1).

ROBIN is intended to support research and validation of rapid workflows. Its outputs require expert interpretation and should not be treated as standalone clinical results.

---

## What ROBIN provides

<div class="grid" markdown>

<div class="robin-feature-card" markdown>
### Real-time processing
ROBIN watches aligned BAM output as sequencing progresses and schedules enabled analyses incrementally.
</div>

<div class="robin-feature-card" markdown>
### Complementary analyses
Methylation classifiers, CNV, MGMT, target/variant analysis and structural-event detection can be combined in one workflow.
</div>

<div class="robin-feature-card" markdown>
### Browser-based monitoring
The web interface brings sample metadata, progress and accumulated analysis results together during a run.
</div>

</div>

---

## Partners and adoption

ROBIN is developed by the **Loose Lab** at the University of Nottingham with collaborators at **Nottingham University Hospitals NHS Trust** and other centres.

- [GitHub repository](https://github.com/LooseLab/ROBIN)
- [Contact form](https://forms.gle/kdX2eiPQPdDUpaBE9)

---

## Information for patients and families

ROBIN is a research tool and is not a substitute for medical advice. General information and support are available from:

- [Cancer Research UK — brain tumours](https://www.cancerresearchuk.org/about-cancer/brain-tumours)
- [The Brain Tumour Charity](https://www.thebraintumourcharity.org/)
- [brainstrust](https://www.brainstrust.org.uk/)
- [Brain Tumour Research](https://www.braintumourresearch.org/)

For concerns about symptoms or treatment, contact an appropriate healthcare professional.
