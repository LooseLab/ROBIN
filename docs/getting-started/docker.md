# Docker (experimental)

!!! abstract "What this page covers"
    Run ROBIN in a container: one application image, bind-mounted BAM and work directories, and the NiceGUI monitor on port **8081**.  
    This is a first-cut packaging of the conda install. Nested tools that themselves call Docker (ClairS-To, MNP-Flex) need extra host-path setup.

!!! warning "Research use only"
    Starting the container with `ROBIN_RESEARCH_ACK=I agree` records the same research-use acknowledgment as typing **I agree** on the host. Do not use this in clinical care.

---

## What is in the image

| Layer | Contents |
|-------|----------|
| Base | `micromamba` (Debian bookworm) |
| Conda | [`robin.yml`](https://github.com/LooseLab/ROBIN/blob/main/robin.yml) — Python 3.12, samtools, bedtools, R/Bioconductor, SnpEff/SnpSift |
| pip | ROBIN from this repository, including git dependencies (Sturgeon, methylartist) |
| Ports | **8081** NiceGUI, **8265** Ray dashboard |

One container is enough for the orchestrator, Ray workers, and the web UI. Separate classifier containers are not required: those jobs already run as Python/R inside ROBIN, or as **sibling** Docker runs (ClairS-To, MNP-Flex) against the host daemon.

---

## Quick start

Compose starts `robin workflow --toml` with [`rcns2_new.workflow-settings.toml`](https://github.com/LooseLab/ROBIN/blob/main/rcns2_new.workflow-settings.toml): the rCNS2 panel, NUH center, the same job list (ITD, CNV, fusion, target, MGMT, classifiers), CNV gene list / penalties, and ITD hotspot mode. The host `work_dir` and reference FASTA are bind-mounted at those same absolute paths.

From a clone with submodules initialized:

```bash
cp docker/compose.env.example .env
# set ROBIN_ADMIN_PASSWORD
# confirm ROBIN_WORK and ROBIN_REFERENCE match the TOML
docker compose up --build
```

Open **http://localhost:8081** and sign in as **admin** with `ROBIN_ADMIN_PASSWORD`.

The TOML has no `path`, so add a watch directory from the GUI after start (same as a native run of this file). Outputs go to `work_dir` on the host (`REF_SAMPLES_NEW` in the current TOML).

Download models (once per image / volume layout):

```bash
docker compose --profile setup run --rm models
# or, after the app is up:
docker compose exec robin robin utils update-models
docker compose exec robin robin utils update-clinvar
```

```bash
docker compose exec robin robin --help
docker compose exec robin robin list-job-types
```

---

## Environment

| Variable | Purpose |
|----------|---------|
| `ROBIN_RESEARCH_ACK` | Must be exactly `I agree` for non-interactive workflow start |
| `ROBIN_ADMIN_PASSWORD` | Creates the first `admin` user when `security.db` is empty |
| `ROBIN_TOML` | Host workflow TOML (default `./rcns2_new.workflow-settings.toml`), mounted at `/config/workflow.toml` |
| `ROBIN_WORK` | Host `work_dir`; must match the TOML and is mounted at that same path |
| `ROBIN_REFERENCE_DIR` | Host directory that contains the FASTA and `.fai`; must cover the TOML `reference` path |
| `GITHUB_TOKEN` | Needed for private model assets |
| `ROBIN_SKIP_MODEL_CHECK` | Set to `1` to silence the entrypoint model warning |
| `LD_LIBRARY_PATH` | Set to `/opt/conda/lib` so GUI sqlite3/ICU use conda's `libstdc++` |
| `XDG_CONFIG_HOME` | GUI auth DB (`security.db`). Compose sets `/home/mambauser/.local` so sqlite is not pointed at the root-owned `~/.config` volume mount |

The GUI password prompt and the terminal **I agree** prompt need a TTY on the host. The env vars above are how a container starts without one.

---

## Nested Docker (ClairS-To, MNP-Flex)

ROBIN already launches some tools with `docker run`. From inside the ROBIN container that only works if:

1. The **host** Docker socket is mounted (`/var/run/docker.sock`).
2. BAM and work paths are bind-mounted at the **same absolute path** on the host and in the ROBIN container. Otherwise the inner container sees the host filesystem and cannot find `/data` or `/work`.

The default compose file already bind-mounts `ROBIN_WORK` and `ROBIN_REFERENCE` at the same absolute host paths used in the TOML. Uncomment the Docker socket volume when you want ClairS-To / MNP-Flex.

Without the socket, SNP calling and MNP-Flex stay unavailable; coverage, target, fusion, ITD, and in-process classifiers still run.

`[minknow] host = "localhost"` in the TOML talks to MinKNOW **inside** the container. To reach MinKNOW on the host, change that host to `host.docker.internal` or run with host networking.

---

## MinKNOW

The container does not include MinKNOW. To talk to a sequencer on the host, use host networking or publish the MinKNOW gRPC port, and install the `minknow` extra in a custom image if you need `robin minknow`. Live watch of BAM folders on a shared mount is the simpler path.

---

## Limits

- Image is large (conda + R + PyTorch). First build takes a long time.
- Models and ClinVar are **not** baked in; download them after build.
- `pywebview` is unused; the UI is the browser.
- This image is Linux/amd64. Apple silicon can run it via emulation but is not a supported deploy path.
- Do not treat this as a production distribution yet.

---

## Files

| Path | Role |
|------|------|
| `Dockerfile` | Application image |
| `docker-compose.yml` | GUI workflow service + optional `models` profile |
| `docker/entrypoint.sh` | Micromamba wrapper + model warning |
| `docker/compose.env.example` | Copy to `.env` (paths aligned with `rcns2_new.workflow-settings.toml`) |
