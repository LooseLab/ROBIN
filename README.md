# ![ROBIN_logo_small.png](src/robin/gui/images/ROBIN_logo_small.png) R.O.B.I.N

**Rapid nanopOre Brain intraoperatIve classificatioN**

> **Research use only.** ROBIN is under active development and validation. Analysis outputs require expert interpretation and must not be treated as standalone clinical results.

ROBIN is a real-time analysis and web-monitoring platform for Oxford Nanopore sequencing of CNS tumours. It processes aligned BAM files as they are produced and can combine methylation classification, copy-number analysis, MGMT analysis, targeted analysis and structural-event detection in one workflow.

ROBIN incorporates **LITTLE JOHN** (Lightweight Infrastructure for Task Tracking and Logging with Extensible Job Orchestration for High-throughput aNalysis) for workflow orchestration and scaling.

The ROBIN approach is described in *Neuro-Oncology*: [Rapid nanopore brain intraoperative classification](https://academic.oup.com/neuro-oncology/article/27/8/2035/8139084).

## Documentation

The MkDocs site is the canonical user and developer documentation:

- [Installation](docs/getting-started/installation.md)
- [Quickstart](docs/getting-started/quickstart.md)
- [MinKNOW configuration](docs/getting-started/minknow-configuration.md)
- [Using ROBIN](docs/using-robin/index.md)
- [Analysis pipelines](docs/analyses/index.md)
- [Command-line reference](docs/cli/index.md)
- [Developer architecture](docs/architecture/overview.md)

## Quick installation

ROBIN requires Python 3.12 and is intended to be installed in a fresh conda environment.

```bash
git clone --recursive https://github.com/LooseLab/ROBIN.git
cd ROBIN

conda env create -f robin.yml
conda activate robin

pip install -e .

robin utils update-models
robin utils update-clinvar
```

Then verify the installation:

```bash
robin --help
robin list-job-types
```

For environment troubleshooting and optional extras, see the [installation guide](docs/getting-started/installation.md).

## Important sequencing requirements

ROBIN's supported real-time workflow expects aligned Oxford Nanopore BAM files produced by MinKNOW.

In particular:

- use HAC basecalling or better;
- enable 5mC/5hmC calling in CpG context when methylation analyses are required;
- perform alignment in MinKNOW using the same reference supplied to ROBIN;
- configure BAM rollover by **read count**, not time;
- keep every input BAM at **50,000 reads or fewer**.

See [MinKNOW configuration](docs/getting-started/minknow-configuration.md) for the complete setup.

## Minimal workflow example

```bash
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/robin_ref/hg38.fa \
  --center Sherwood \
  --target-panel rCNS2
```

A recommended first-run sequence, including reference staging and opening the GUI, is in the [Quickstart](docs/getting-started/quickstart.md). Full command options are documented in the [CLI reference](docs/cli/index.md).

## System requirements

A typical production-scale installation should have approximately **64 GB RAM or more**. CPU and GPU requirements depend on the Oxford Nanopore sequencing and basecalling configuration in use.

Some optional downstream analyses require additional software such as Docker. ROBIN's core scientific, bioinformatics and R dependencies are defined in `robin.yml` and `pyproject.toml`.

## Repository layout

- `src/robin/` — application, workflow engine, analyses and GUI
- `docs/` — MkDocs user and developer documentation
- `robin.yml` — primary conda environment
- `pyproject.toml` — Python package configuration
- `src/robin/resources/` — packaged panels, manifests and other resources

## Reporting problems

Please use [GitHub Issues](https://github.com/LooseLab/ROBIN/issues) for reproducible software problems and feature requests.

When reporting a workflow problem, include the ROBIN version/commit, operating system, command used, relevant log output and enough information about the input configuration to reproduce the issue without sharing identifiable patient data.

## License

ROBIN is distributed under the **CC BY-NC 4.0** license. See [LICENSE](LICENSE).

Third-party tools, models and resources retain their own licences and terms of use.

## Acknowledgments

ROBIN builds on work from the nanopore and computational pathology communities, including [Sturgeon](https://github.com/marcpaga/sturgeon), [Rapid-CNS2](https://link.springer.com/article/10.1007/s00401-022-02415-6), [Readfish](https://github.com/LooseLab/readfish), `cnv_from_bam` and `methylartist`.

ROBIN is developed by the **Loose Lab, University of Nottingham**, with collaborators at **Nottingham University Hospitals NHS Trust** and other centres.
