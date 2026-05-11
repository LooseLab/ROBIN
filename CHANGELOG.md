# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project (almost) adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- **Job batching (workflow_ray.py & workflow_simple.py):** Reduced fragmentation of analysis jobs into many small batches when BAMs arrive slower than the previous fixed 2 s batcher timeout.
  - `SampleJobBatcher` now uses an adaptive timeout: an "idle" timeout when no jobs of that type are running on a worker (default **2 s**, keeps the first batch reactive) and a longer "busy" timeout while at least one job of that type is in flight (default **30 s**, lets more files accumulate into the next batch).
  - When a worker for a given job type frees up, the batcher is force-flushed for that type so the worker has something to pick up immediately rather than waiting out the busy timeout.
  - The dispatcher now coalesces compatible waiting batched jobs at dispatch time. When a batched job is about to be dispatched (or a worker pulls one from a per-type queue), other waiting batched jobs for the same `(sample_id, job_type)` are merged into it up to `max_batch_size`. Jobs without `_batched_job` metadata or any context tagged `force_individual_batch` (e.g. large BAMs) are never coalesced.
  - New env vars to tune the per-type adaptive timeouts:
    - `ROBIN_BATCH_TIMEOUT_IDLE_S` — global idle timeout (seconds)
    - `ROBIN_BATCH_TIMEOUT_BUSY_S` — global busy timeout (seconds)
    - `ROBIN_BATCH_TIMEOUT_IDLE_S_<TYPE>` / `ROBIN_BATCH_TIMEOUT_BUSY_S_<TYPE>` — per-type override (e.g. `ROBIN_BATCH_TIMEOUT_BUSY_S_CNV=45`)
  - Coordinator `stats()` (Ray) and `get_stats()` (simple) now expose `coalesce_events`, `coalesce_jobs_absorbed`, and `coalesce_contexts_merged` for observability.

## [0.5.2] - 2026-05-05

### Changed
- Improved workflow throughput for very large datasets by refining queue backpressure and job-admission behavior so the pipeline sustains progress more reliably under heavy load.

## [0.5.1] - 2026-04-30

### Changed
- Removed **PanCan** as a built-in target panel across CLI, analysis, GUI panel mappings, and documentation.
- Removed packaged PanCan panel resources: `PanCan_panel_name_uniq.bed` and `PanCan_panel_source.bed`.
- Kept **PanNanoDX** classifier support unchanged (`pannanodx` workflows and model references remain available).
- Bumped project/app version metadata from **`0.5`** to **`0.5.1`**.

## [0.5] - 2026-04-04

### Changed
- **Conda environment:** a single **`robin.yml`** (env name **`robin`**) replaces separate `robin_0_8.yml`, `robin_py12.yml`, `robin_osx.yml`, and the old Python 3.9 `robin.yml` layout. Create with `conda env create -f robin.yml` and `conda activate robin`.
- **Linux `libstdc++`:** optional `conda_env_hooks/` activate scripts and `scripts/setup_conda_ld_path.sh` were removed; use **`export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${LD_LIBRARY_PATH}"`** (see README *Common issues*) and, on Linux, **`robin_linux_extras.yml`** or `conda install` for `libstdcxx-ng>=13` / `libgcc-ng>=13` (plain `conda env create` does not honour `# [linux]` selectors in `robin.yml`).

### Important
- **Main branch merge preparation:** This release line captures the major merge of the `little_john` integration work back into the main repository branch.
- **Version update:** Project version is now **`0.5`**.
- **Environment reset required:** **Create a brand-new conda environment for `0.5`. Do not reuse any previous ROBIN/little_john conda environment.**

### Added
- **Watched folders (smart add / Option B):** When adding a folder to watch, ROBIN scans for mixed “already analysed” vs “new” samples and automatically adds only the subfolders that contain exclusively un-analysed samples. Previously-analysed samples are skipped and reported.
- **Watched folders (GUI feedback):** If adding a folder results in some samples/subfolders being skipped, the GUI now shows a warning notification rather than a plain success.
- **Job type filtering:** Added job-type filtering to improve workflow/GUI navigation for large runs.
- **SNP calling global controls:** Added global GUI controls/buttons for SNP calling to simplify starting/monitoring variant calling.
- **MNP-Flex global controls:** Added global GUI controls/buttons for MNP-Flex workflows.
- **Process large BAMs option:** When `ROBIN_PROCESS_LARGE_BAMS` is set (e.g. `1`, `true`, `yes`, `on`), BAMs with more than 50,000 reads are no longer skipped. They are processed with per-file batching so each large BAM is sent to workers as a single-file batch (no batching with other files). Default remains to skip BAMs over 50k reads.
- **Launch warnings for large-BAM mode:** When `ROBIN_PROCESS_LARGE_BAMS` is enabled, the GUI shows a warning in the disclaimer dialog (if shown) and a one-time notification on launch; the CLI prints a warning after the disclaimer when starting the workflow. Users are advised not to use this option alongside live runs.
- **Fusion GUI delay analysis:** `FUSION_GUI_DELAY_ANALYSIS.md` documenting the cause of the ~29s delay between Summary and Fusion tab updates (30s refresh timer) and the mitigation (immediate refresh on section build).
- **Fusion GUI timing instrumentation:** `refresh_fusion()` and `_load_fusion_data()` in the Fusion component now log duration (load time, UI update time, total) to help diagnose slow fusion processing (e.g. `refresh_fusion() completed in X.XXs (load=Y.YYs, ui=Z.ZZs)`).
- **GUI authentication:** All GUI routes require login. Unauthenticated users are redirected to a themed `/login` page. Logout is available from the hamburger menu only (styled like the Quit button).
- **Password stored as Argon2 hash** in `~/.config/robin/gui_password_hash` (or `%APPDATA%/robin/gui_password_hash` on Windows). No plaintext or env-var password.
- **Password set at GUI startup:** On first run, prompt twice (set + confirm) with hidden input; on later runs, prompt once to verify. Skipped when not a TTY (e.g. headless); users then log in via the web login page.
- **`robin password set`** CLI to set or replace the GUI password. Asks to replace if a hash already exists; uses Rich for prompts when available. Password input is never echoed.
- **Session invalidation on server restart:** A per-run auth generation token ensures persisted user storage from a previous run (or before a password reset) is not trusted; users must log in again after restart or password change.
- Dedicated "Manage watched folders" page at `/watched_folders` (replaces modal); accessible via link from the welcome page.
- Loading modal when adding a watched folder: shows "Adding folder..." with spinner and closes automatically on success or error.
- Path overlap validation in `add_watch_path`: rejects paths that equal, contain, or are inside the work directory to prevent watching already-analysed output.
- `parquet_filter.txt` as preferred CpG filter for parquet creation (falls back to `sturg_nanodx_cpgs_0125.bed.gz` when absent).
- `scripts/compare_mgmt_parquet_methylation.py` to compare MGMT methylation counts between bed and parquet for validation.
- `METHYLATION_EXTRACTION_COMPARISON.md` documenting differences between MNP-Flex (parquet) and MGMT (bedmethyl) methylation extraction, including coordinate handling and BAM set divergence.
- SNP analysis pipeline via Clair3 with snpEff and SnpSift annotation.
- SNP display data generation and GUI table with IGV navigation. Current extensive logging to the command line to track progress.
- IGV-ready BAM creation for consistent genome browser loading.
- Workflow support for `igv_bam` and `snp_analysis` job types.
- SNP analysis queueing from the GUI with concurrency control.
- MNP-FLEX upload support for bedMethyl outputs for beta users with Epignostix credentials (https://epignostix.com/). Supply credentials via `MNPFLEX_USER` and `MNPFLEX_PASS` environment variables.
- `robin utils mgmt` command to summarize MGMT CpG site methylation counts from `mgmt_sorted.bam` outputs.

### Changed
- **Design refresh:** Updated GUI styling/layout (design refresh).
- **Classification scoring:** Improved/fixed classification scoring and aggregation across classifiers (Sturgeon/NanoDX/PanNanoDX/RandomForest).
- **Background work / scheduling:** Improved background processing behavior (responsiveness/robustness).
- **Watched folders (performance):** Folder scanning during “Add folder to watch” now probes only a small number of BAMs per directory (sampling) to keep the GUI responsive for large folders.
- **Batching for large BAMs:** When `ROBIN_PROCESS_LARGE_BAMS` is set, the preprocessor marks large BAMs (>50k reads) with `force_individual_batch`. The Ray and simple workflow batchers now emit single-file batches for such jobs and batch remaining jobs as before, so large BAMs are always processed one per worker.
- **Report layout and typography:**
  - **Classification:** Methylation classification plots now use seaborn for improved styling; plots show top 3 classes with larger fonts and thicker lines. Detailed classification tables use a 2-column layout (Sturgeon|NanoDX, PanNanoDX|Random Forest) with compact styling and top 10 predictions.
  - **CNV:** Summary and event tables use compact styling with unified 9pt font. Chromosome Arm Events and CNV Events Containing Genes tables are displayed side-by-side in two columns. Arm/Region columns widened to prevent text wrapping (e.g. "p-arm"). Increased cell padding and line spacing for readability.
  - **Coverage:** Removed overlapping outlier labels from the "Coverage Distribution by Chromosome" box plot. Added a "Coverage outliers (potential gains/losses)" table below the plot with Chromosome, Gene, Coverage, and Type (Gain/Loss) columns.
  - **Base:** `create_table()` now supports `compact`, `font_size`, and improved leading for consistent table typography across sections.
- Parquet filter loading: `.txt` files (e.g. `parquet_filter.txt`) are treated as 1-based Illumina coordinates and converted to 0-based for BED/PyRanges; fixes off-by-one in CpG site filtering.
- Parquet filter cache: distinct cache suffix (`_1based`) for `.txt` filters to avoid stale data after coordinate conversion change.
- Parquet filter format: support for header row and whitespace-separated columns in `parquet_filter.txt`.
- **Note:** After the parquet filter coordinate fix, existing parquet files should be rebuilt (delete parquet and re-run bed_conversion) to obtain correct methylation counts.
- IGV viewer initialization and BAM loading prioritization.
- Target analysis now always emits `targets_exceeding_threshold.bed` for SNP analysis.
- MGMT methylation classification uses per-read probabilities (fixes aggregated max bleed).
- Version metadata aligned to **`0.5`** across packaging and app entry points.
- `robin utils mgmt` now reports methylation percent from counts and normalizes bedmethyl fraction inputs.
- MNP-Flex results panel now hides empty fields until results are available.
- `robin utils mgmt` supports recursive search and file output for TSV summaries.
- Command-line progress output supports Rich styling; updating the installation is required to enable Rich functionality.
- **Fusion GUI:** Fusion tab now refreshes immediately when the Fusion section is built (in addition to the 0.5s deferred and 30s periodic timers), so fusion data appears in sync with the Summary without waiting for the next timer tick.
- **Ray per-job memory (workflow_ray.py):** Configurable memory per job type — 1 GiB default for most jobs (preprocessing, bed_conversion, mgmt, cnv, target, classification, igv_bam, target_bam_finalize), 4 GiB for fusion, 8 GiB for snp_analysis (Clair3/variant calling). Reduces memory pressure on shared nodes while keeping SNP calling at 8 GiB.
- **CNV plots (GUI):** Bin width selector added to the CNV plots so users can choose the bin width used for visualization (independent of the analysis bin width).
- **CNV plots (GUI):** Chromosome filtering now matches the report: only chr0–chr22, chrX, and chrY are shown in the chromosome selector and in the "All" view. chrM and other contigs (e.g. unplaced, alt) are excluded.
- Bump version to 0.5 (major backend upgrade to Python 3.12).
- **Quit button (menu):** Shown only when the app is opened from localhost or 127.0.0.1 (same device as the server). Hidden when visiting via any other host (e.g. remote or LAN URL) so only someone at the sequencing machine can quit the app.
- **Python 3.12:** Minimum supported Python version is now 3.12.

### Fixed
- **GUI locking / concurrency:** Fixes to GUI locking to prevent UI/workflow concurrency issues.
- **Timer/refresh reliability:** Fixes to timer-driven refresh/update behavior.
- **Multi-BAM samples:** Removed per-BAM “already analysed” skipping during preprocessing to avoid accidentally skipping valid BAMs from the same sample/run.
- **Ray dedup for SNP analysis:** Added SNP calling to the dedup list to prevent duplicate SNP runs for a sample.

### Dependencies
- Added argon2-cffi for GUI password hashing.
- Added seaborn for methylation classification plot styling.
- Added Rich for enhanced CLI progress and styled output (requires updating installation).
- Updated bundled ClinVar resource (ClinVar VCF + index).

## [0.0.2] - 2024-12-19

### Added
- Model setup scripts and asset management system
- Performance optimizations for reporting system
- Enhanced error handling and loading fixes

### Changed
- Migrated from Git LFS to release asset system for model files
- Updated README with correct repository URL and installation instructions
- Improved performance across multiple analysis modules

### Fixed
- Fusion detection analysis fixes
- Tag error handling improvements
- BED file parsing to allow '.' in column 5
- Minimum read count validation (ensures minimum of 3 reads)
- Various error handling improvements

### Performance
- Speed improvements in reporting system
- Multiple performance tweaks across analysis modules
- Optimized loading and processing workflows

## [0.0.1] - 2024-12-14

### Added
- Initial release of robin bioinformatics workflow engine
- BAM file processing and analysis capabilities
- Multi-queue workflow system
- Built-in analysis modules (MGMT, CNV, fusion detection, etc.)
- GUI workflow monitor using NiceGUI
- Ray distributed computing support
- Comprehensive documentation and installation guides
