# MARLIN reference assets

Bundled reference files for the MARLIN leukemia methylation classifier
(https://github.com/hovestadt/MARLIN, MIT).

| File | Purpose |
|------|---------|
| `marlin_v1.features.RData` | Ordered 357,340 CpG probe IDs |
| `marlin_v1.class_annotations.xlsx` | 42 class labels / metadata |
| `marlin_v1.probes_{hg19,hg38,t2t}.bed.gz` | Probe genomic coordinates |

The trained Keras weights (`marlin_v1.model.hdf5`, ~1.1 GiB) are **not**
vendored. ROBIN downloads them on first MARLIN job from Zenodo record
[15565404](https://zenodo.org/records/15565404) into `~/.cache/robin/marlin/`.

Overrides:

- `ROBIN_MARLIN_MODEL_PATH` — explicit model file
- `ROBIN_MARLIN_CACHE_DIR` — cache directory for auto-download
- `ROBIN_MARLIN_MODEL_URL` — alternate download URL

Install via `pip install 'robin[marlin]'` (TensorFlow ≥2.16 + `tf-keras` on Python 3.12+).
Inference uses the Keras 2-compatible `tf_keras` loader because the MARLIN
checkpoint is a legacy HDF5 saved with Keras 2.
