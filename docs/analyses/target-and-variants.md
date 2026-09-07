# Target coverage and variants

Target analysis operates on the active ROBIN target panel. It provides targeted coverage information and prepares data used by downstream variant-analysis paths.

The main implementation is in `src/robin/analysis/target_analysis.py`.

## Target panels

A target panel defines genomic intervals/genes of interest. ROBIN ships with supported built-in panels and can register custom BED-based panels through the CLI.

Use:

```bash
robin list-panels
```

to inspect available panels. See [Panel commands](../cli/panels.md) for adding and removing custom panels.

## Coverage analysis

For each sample, target analysis accumulates aligned-read evidence over panel regions and generates target-level coverage summaries. These outputs support both direct interpretation in the GUI and downstream steps that need to know which regions have sufficient sequencing evidence.

Because sequencing is continuous, coverage should be viewed as a time-dependent quantity: a target below threshold early in a run may exceed the threshold later.

## Variant-analysis path

`target_analysis.py` also contains support for downstream SNP analysis. The current implementation includes:

- preparation of target BAM data;
- Clair3-based variant calling;
- snpEff annotation;
- SnpSift annotation against ClinVar where the required database/index is available;
- chromosome-name normalisation for ClinVar annotation;
- generation of display-ready variant data.

The variant path has additional dependencies beyond target coverage itself.

## Docker requirement

The Clair3 path uses Docker. ROBIN checks whether the Docker Python package is installed and whether the Docker daemon is reachable before attempting this analysis.

If Docker is unavailable, target coverage can still be useful, but the Docker-dependent variant step will not run successfully.

## ClinVar resources

ROBIN expects a bgzipped ClinVar VCF and tabix index for SnpSift annotation. Assets can be refreshed with:

```bash
robin utils update-clinvar
```

The implementation attempts to ensure that the tabix index exists before annotation.

## Reference and chromosome naming

Variant annotation frequently combines resources that use different contig naming conventions. For example, one VCF may use `chr1` while another uses `1`.

ROBIN includes a normalisation/restoration step around ClinVar annotation to handle this mismatch, but all primary sequencing/reference resources should still be based on the same genome assembly.

## Interpretation

Coverage and variant calls answer different questions:

- coverage indicates how much sequencing evidence exists over a target;
- a variant caller evaluates sequence differences supported by those reads;
- annotation adds external information to called variants.

A region having high coverage does not imply that a variant is present, and a variant observed at low coverage requires particular caution.

## Troubleshooting

### Target panel not recognised

Run `robin list-panels` and verify the panel name supplied to `--target-panel`. For custom panels, confirm registration succeeded and the BED coordinates use the expected reference.

### Variant calling does not start

Check Docker first:

```bash
docker info
```

Also verify that the current user has permission to access the Docker daemon and that the required image can be obtained.

### ClinVar annotation missing

Refresh ClinVar resources with `robin utils update-clinvar` and check that both the compressed VCF and `.tbi` index are present in the ROBIN resources location.

### Coverage differs from expectation

Confirm the BAM/reference/panel coordinate systems match, then inspect individual alignments in the affected target. Differences can also arise from mapping quality, supplementary/secondary alignments, or incomplete accumulated sequencing depth.