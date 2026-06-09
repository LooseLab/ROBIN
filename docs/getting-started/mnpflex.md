# MNP-Flex setup

ROBIN can submit methylation data to **MNP-Flex** and display the returned
classification, quality-control, and MGMT results in the web interface.
MNP-Flex is a commercial service provided by
[Heidelberg Epignostix GmbH](https://epignostix.com/).

## Before you begin

You need:

- An agreement with Epignostix covering use of MNP-Flex.
- An Epignostix username and password with access to the MNP-Flex service.
- A ROBIN workflow that produces the per-sample methylation parquet file.
- Network access from the machine running ROBIN to the Epignostix service.

ROBIN sends sample-derived methylation data to an external service. Confirm
that this is permitted by your local information-governance and data-transfer
policies before enabling the integration.

## Configure credentials

Set the credentials in the environment of the **server process running
ROBIN**, before starting `robin workflow`:

```bash
export MNPFLEX_USERNAME="your-epignostix-username"
export MNPFLEX_PASSWORD="your-epignostix-password"
```

Then start ROBIN from the same shell:

```bash
robin workflow <data_folder> \
  --work-dir <output_folder> \
  --reference <reference.fa> \
  --center <center_id>
```

Both variables are required. If either is absent, ROBIN hides the MNP-Flex
controls and does not submit samples.

Do not put credentials in the repository, command history, workflow output,
screenshots, or other files that may be shared. For a persistent deployment,
configure these variables through the service manager or secret-management
system used to launch ROBIN.

!!! note "Legacy variable names"

    Older ROBIN code and release notes may refer to `MNPFLEX_USER` and
    `MNPFLEX_PASS`. Current GUI and batch integration uses
    `MNPFLEX_USERNAME` and `MNPFLEX_PASSWORD`; use the current names for new
    deployments.

## Optional settings

The defaults match the standard Epignostix integration. Change these only when
Epignostix supplies different values:

| Variable | Default | Purpose |
|----------|---------|---------|
| `MNPFLEX_BASE_URL` | `https://app.epignostix.com` | Epignostix API base URL |
| `MNPFLEX_WORKFLOW_ID` | `18` | MNP-Flex workflow identifier |
| `MNPFLEX_CLIENT_ID` | `ROBIN` | OAuth client identifier |
| `MNPFLEX_CLIENT_SECRET` | `SECRET` | OAuth client secret |
| `MNPFLEX_SCOPE` | empty | Optional OAuth scope |

For example:

```bash
export MNPFLEX_WORKFLOW_ID="18"
export MNPFLEX_BASE_URL="https://app.epignostix.com"
```

## Confirm that it is enabled

After restarting ROBIN:

1. Open a sample page.
2. Look for the **MNP-Flex results** section.
3. Use **Generate MNP-Flex subset BED** to confirm that the sample methylation
   parquet can be prepared.
4. Use **Run MNP-Flex analysis** to submit that sample.

The samples overview also shows **mnpflex run all** when credentials are
available. This runs MNP-Flex sequentially for eligible samples that do not
already have results.

ROBIN may automatically submit a completed sample when its analysis jobs have
finished and no MNP-Flex result is present.

## Inputs and outputs

ROBIN reconstructs bedMethyl data from the sample parquet, prepares the subset
required by MNP-Flex, uploads it, and retrieves the result bundle.

Files are written within the sample output directory:

```text
<sample-id>/
├── <sample-id>.mnpflex.bed
├── <sample-id>.MNPFlex.subset.bed
└── mnpflex_results_<sample-id>/
    └── bundle_summary.json
```

Additional plots and result files may be present in the results directory.
These files feed both the sample page and ROBIN-generated reports.

## Troubleshooting

**The MNP-Flex section or bulk button is missing**

- Confirm both `MNPFLEX_USERNAME` and `MNPFLEX_PASSWORD` are set for the
  process that launches ROBIN.
- Restart ROBIN after changing its environment.
- If ROBIN is launched by a service, setting variables in an interactive shell
  does not change the service environment.

**The sample cannot generate a subset BED**

- Confirm the sample output directory contains `<sample-id>.parquet` or
  another sample parquet file.
- Confirm methylation preprocessing completed successfully.

**Authentication or upload fails**

- Check the credentials with Epignostix.
- Confirm the machine can reach `MNPFLEX_BASE_URL`.
- Confirm any custom workflow ID, client ID, client secret, or scope with
  Epignostix.
- Review the ROBIN log for messages containing `MNPFlex`.

**Results do not appear immediately**

- MNP-Flex processing is asynchronous and can take time.
- Use the section's status and refresh controls, and check the ROBIN log for
  upload, polling, or result-download errors.
