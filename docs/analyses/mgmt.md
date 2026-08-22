# MGMT promoter methylation

ROBIN includes a focused analysis of methylation in the **MGMT** promoter region. The workflow implementation is in `src/robin/analysis/mgmt_analysis.py`.

## Region analysed

The current implementation targets:

```text
chr10:129466536-129467536
```

Coordinates must be interpreted against the reference assembly expected by ROBIN. Input BAMs should use compatible chromosome names and coordinates.

## Analysis path

The MGMT module combines locus extraction, modified-base information, model prediction, and visualisation.

At a high level it:

1. checks whether the incoming BAM contains reads overlapping the MGMT locus;
2. extracts/maintains the relevant alignments;
3. safely sorts and indexes the locus-specific BAM data;
4. evaluates methylation at the required sites;
5. applies the configured R-based prediction model;
6. generates sample-level result data and visualisation for the GUI/reporting layer.

The implementation includes explicit checks that sorted BAMs are complete and readable before indexing, reducing the risk that a partially written file is propagated into downstream analysis.

## Input requirements

MGMT analysis requires:

- an aligned BAM using the expected reference coordinates;
- reads covering the MGMT promoter region;
- modified-base calls suitable for CpG methylation analysis;
- the ROBIN model/R resources installed with the environment.

If there are no reads at the locus, the absence of an MGMT prediction should not be interpreted as a biological negative result.

## Real-time behaviour

As more BAM chunks arrive, additional locus-spanning reads can contribute to the sample-level analysis. Confidence in the methylation estimate therefore depends on accumulated informative coverage.

Early results should be interpreted cautiously where only a small number of reads/CpG observations are available.

## Interpretation

The MGMT result is a model-based methylation prediction. Review it together with the underlying methylation evidence and other sample information.

Important considerations include:

- number of informative reads;
- consistency of methylation across relevant CpGs;
- sequencing/basecalling quality;
- whether the expected modified-base model was used;
- whether the reference coordinates match the BAM.

## Troubleshooting

### No reads at the MGMT locus

Confirm that the BAM is aligned, contains `chr10`, and uses the same reference build expected by ROBIN. You can independently inspect the locus with `samtools view` or IGV if necessary.

### BAM sorting/indexing errors

The MGMT implementation validates intermediate BAMs before indexing. Errors at this stage can indicate a truncated/incomplete input file or filesystem problem. Check the job log and verify the source BAM with `samtools quickcheck`.

### No methylation prediction despite coverage

Check that modified-base tags are present and that CpG-context 5mC/5hmC calling was enabled during basecalling. Also confirm that the required R/model resources are available in the active ROBIN environment.