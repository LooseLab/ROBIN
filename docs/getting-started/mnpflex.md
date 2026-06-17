# MNP-Flex setup

ROBIN can run **MNP-Flex** locally via Docker or submit methylation data to the
**Epignostix API**, then display classification, quality-control, and MGMT
results in the web interface.

## Before you begin

You need:

- A ROBIN workflow that produces the per-sample methylation parquet file.
- Either a loaded MNP-Flex Docker image **or** Epignostix API credentials.
- For the API path: network access from the machine running ROBIN to Epignostix.

For Docker, confirm that running MNP-Flex locally is permitted by your local
information-governance policies. For the API, confirm that sending
sample-derived methylation data to an external service is permitted.

!!! warning "Wait for sufficient sequencing data"

    MNP-Flex is recommended only after at least **12 hours of sequencing data**
    have been generated for the sample. Running it earlier may provide less
    reliable or less representative results.

## Choose a backend

Set `MNPFLEX_BACKEND` in the environment of the **server process running
ROBIN**, before starting `robin workflow`:

| Backend | Required settings |
|---------|-------------------|
| `docker` | `MNPFLEX_DOCKER_IMAGE` |
| `api` | `MNPFLEX_USERNAME`, `MNPFLEX_PASSWORD` |
| `disabled` | Hide analysis runs (BED generation remains available) |

If `MNPFLEX_BACKEND` is not set, ROBIN keeps legacy behaviour: API credentials
alone enable the Epignostix backend.

### Local Docker

```bash
export MNPFLEX_BACKEND=docker
export MNPFLEX_DOCKER_IMAGE=mnpflex-synnovis:1.0.0

# Load the image once (example Synnovis bundle)
docker load < mnpflex-synnovis-1.0.0.tar.gz
docker images mnpflex-synnovis
```

Use the exact image name and tag shown by `docker images`. Different sites may
use different registry paths or versions, for example
`registry.example.com/mnpflex:2.1.0`.

Optional Docker settings:

| Variable | Default | Purpose |
|----------|---------|---------|
| `MNPFLEX_DOCKER_TIMEOUT` | `3600` | Max seconds per Docker run |
| `MNPFLEX_DOCKER_BINARY` | `docker` | Container CLI (`podman` also works) |
| `MNPFLEX_DOCKER_EXTRA_ARGS` | empty | Extra flags passed to `docker run` |
| `MNPFLEX_DOCKER_INPUT` | `full` | `full` uses `<sample>.mnpflex.bed`; `subset` uses the panel subset BED |

### Epignostix API

```bash
export MNPFLEX_BACKEND=api
export MNPFLEX_USERNAME="your-epignostix-username"
export MNPFLEX_PASSWORD="your-epignostix-password"
```

Optional API settings:

| Variable | Default | Purpose |
|----------|---------|---------|
| `MNPFLEX_BASE_URL` | `https://app.epignostix.com` | Epignostix API base URL |
| `MNPFLEX_WORKFLOW_ID` | `18` | MNP-Flex workflow identifier |
| `MNPFLEX_CLIENT_ID` | `ROBIN` | OAuth client identifier |
| `MNPFLEX_CLIENT_SECRET` | `SECRET` | OAuth client secret |
| `MNPFLEX_SCOPE` | empty | Optional OAuth scope |

!!! note "Legacy variable names"

    Older ROBIN code and release notes may refer to `MNPFLEX_USER` and
    `MNPFLEX_PASS`. Current GUI and batch integration uses
    `MNPFLEX_USERNAME` and `MNPFLEX_PASSWORD`; use the current names for new
    deployments.

Then start ROBIN from the same shell:

```bash
robin workflow <data_folder> \
  --work-dir <output_folder> \
  --reference <reference.fa> \
  --center <center_id>
```

Do not put credentials in the repository, command history, workflow output,
screenshots, or other files that may be shared.

## Confirm that it is enabled

After restarting ROBIN:

1. Open a sample page.
2. Look for the **MNP-Flex results** section and the active backend label.
3. Use **Generate MNP-Flex subset BED** to confirm that the sample methylation
   parquet can be prepared.
4. Once at least **12 hours of sequencing data** have been generated, use
   **Run MNP-Flex analysis**.

The samples overview shows **mnpflex run all** when Docker or API analysis is
configured. ROBIN may automatically submit a completed sample when its analysis
jobs have finished and no MNP-Flex result is present.

## Inputs and outputs

ROBIN reconstructs bedMethyl data from the sample parquet, prepares the BED
input required by the selected backend, runs MNP-Flex, and writes a normalized
result bundle for the GUI and reports.

Files are written within the sample output directory:

```text
<sample-id>/
├── <sample-id>.mnpflex.bed
├── <sample-id>.MNPFlex.subset.bed
└── mnpflex_results_<sample-id>/
    ├── bundle_summary.json
    ├── qc_coverage_plot.png
    ├── qc_methylation_density_plot.png
    ├── mgmt_region_plot.png
    └── docker_raw/              # Docker path only: original CSV/PNG outputs
```

## Troubleshooting

**The MNP-Flex section shows a configuration error**

- For Docker: set `MNPFLEX_BACKEND=docker` and `MNPFLEX_DOCKER_IMAGE`.
- For API: set `MNPFLEX_BACKEND=api` with username and password.
- Restart ROBIN after changing its environment.

**The bulk button is missing**

- Confirm Docker or API analysis is configured and valid.
- Restart ROBIN after changing its environment.

**Docker run fails**

- Confirm `docker info` works for the ROBIN process user.
- Confirm `docker image inspect $MNPFLEX_DOCKER_IMAGE` succeeds.
- Review the ROBIN log for `[MNPFlex] Running Docker`.

**Authentication or API upload fails**

- Check the credentials with Epignostix.
- Confirm the machine can reach `MNPFLEX_BASE_URL`.
- Review the ROBIN log for messages containing `MNPFlex`.

**The sample cannot generate a subset BED**

- Confirm the sample output directory contains `<sample-id>.parquet` or
  another sample parquet file.
- Confirm methylation preprocessing completed successfully.

**Results do not appear immediately**

- API processing is asynchronous and can take time.
- Docker runs are synchronous but may take several minutes on large samples.
- Check the ROBIN log for upload, Docker, or adaptation errors.
