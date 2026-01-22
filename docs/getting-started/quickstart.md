# ROBIN Quickstart

This guide covers quickstart instructions for both versions of ROBIN.

- [Robin & Little John](#robin--little-john)
- [Original Robin](#original-robin)

For installation details, see the [Installation Guide](installation.md).

# Robin & Little John

## Usage

ROBIN expects to analyse BAM files generated during sequencing by an Oxford Nanopore Technologies sequencer. ROBIN presumes real time HAC basecalling (SUP is not required). ROBIN presumes data have been called with 5hmC 5mC methylation calling in MinKNOW. ROBIN presumes real time alignment is running in MinKNOW - ROBIN does not realign your reads.

ROBIN assumes that BAM files are being output in small batches. We recommend setting file output to one bam for every 10,000 to 50,000 reads. We do not support real time processing of BAM files in 1 hour chunks (the default output).

ROBIN does not use pod5 data or fastq data from the sequencer - you can deselect these options if you wish.

## Important

On platforms with 64Gb of RAM or less we recommend restarting your device prior to a run. As an example, if you are running a p2i and start a run on position A and later on position B we recommend you restart at the end of the run on position B (once base calling is complete). You can simply restart dorado if you know how to do this. We find dorado holds on to memory for an indefinite period of time and this can cause problems.

## Run a Workflow

The primary command for running robin workflows is:

```bash
robin workflow <data_folder> --work-dir <output_folder> -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest --reference ~/references/hg38_simple.fa --center <center_id>
```

### Command Breakdown

- `robin workflow`: The main workflow command
- `<data_folder>`: Directory containing your BAM files
- `--work-dir <output_folder>`: Directory where results will be saved
- `-w`: Workflow specification (comma-separated list of analysis types)
- `--reference`: Path to reference genome (required for some analyses)
- `--center <center_id>`: Center ID running the analysis (e.g., `Sherwood`, `Auckland`, `New York`)
- `--target-panel`: The specific panel that is being applied

## Example Usage

```bash
# Basic workflow with all analysis types
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center Sherwood \
  --target-panel rCNS2

# Simplified workflow with just a few analyses
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center Auckland \
  --target-panel PanCan

# With verbose output and custom logging
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,cnv,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center New_York \
  --target-panel rCNS2 \
  --verbose \
  --log-level INFO
```

## Known Issues

- Currently SNP calling is not enabled in this version of ROBIN. It will be re-enabled in the near future.
- CNV change inference is based on extensive heuristics - every call should be checked by visual inspection.
- If you ctrl-c to end ROBIN it will do its best to clean up and stop gracefully but may fail.
- CSV data export is in development but is not currently available - it will be enabled in the near future.
- Many other unknown issues - please open an issue and we will resolve where possible.

## Available Commands

### list-job-types

List all available job types organized by queue category:

```bash
robin list-job-types
```

Available Job Types:

- Preprocessing Queue: preprocessing
- Bed Conversion Queue: bed_conversion
- Analysis Queue: mgmt, cnv, target, fusion
- Classification Queue: sturgeon, nanodx, pannanodx
- Slow Queue: random_forest

### workflow

Run an async workflow on BAM files in a directory:

```bash
robin workflow /path/to/directory --workflow "workflow_plan" [OPTIONS]
```

Required Options:

- `--workflow`, `-w`: Workflow plan (e.g., `mgmt,sturgeon` or `preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon`)
- `--center`: Center ID running the analysis (e.g., `Sherwood`, `Auckland`, `New York`)

Optional Options:

- `--work-dir`, `-d`: Base output directory for analysis results
- `--reference`, `-r`: Path to reference genome (FASTA format)
- `--verbose`, `-v`: Enable verbose output and detailed error traces
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--log-level`: Global log level (DEBUG|INFO|WARNING|ERROR, default: ERROR)
- `--job-log-level`: Set log level for specific job (e.g., `preprocessing:DEBUG`, `mgmt:WARNING`)
- `--deduplicate-jobs`: Job types to deduplicate by sample ID (e.g., `sturgeon`, `mgmt`)
- `--no-progress`: Disable progress bars for file processing
- `--use-ray`/`--no-use-ray`: Enable Ray distributed computing (default: on)
- `--with-gui`/`--no-gui`: Launch NiceGUI workflow monitor (default: on)

## Panel Management

Manage built-in and custom target panels used by analyses like target, cnv, and fusion.

Built-in panels include: rCNS2, AML, PanCan.
Custom panels are stored internally after you add them from a BED file.

### List available panels

```bash
robin list-panels
```

### Add a custom panel from a BED file

```bash
# Add and register a panel (BED must have >=4 columns: chr, start, end, gene_name[s])
robin add-panel /path/to/your_panel.bed MyCustomPanel

# Optional: validate format only, without adding
robin add-panel /path/to/your_panel.bed MyCustomPanel --validate-only
```

Notes:

- Panel names cannot be empty and cannot reuse reserved names: rCNS2, AML, PanCan.
- BED may be 4- or 6-column; if multiple genes are in one region, use comma-separated names.

### Remove a custom panel

```bash
# Will prompt for confirmation
robin remove-panel MyCustomPanel

# Skip confirmation
robin remove-panel MyCustomPanel --force
```

Built-in panels cannot be removed.

### Use a panel in a workflow

```bash
robin workflow /path/to/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion \
  --target-panel MyCustomPanel \
  --center Sherwood
```

## Performance Features

### Enhanced Processing Capabilities

- Batched Processing: All analysis workflows now support batched processing for improved efficiency
- Memory Optimization: Intelligent memory management for large datasets and long-running processes
- Multi-threading: Configurable multi-threaded BAM processing via `LJ_BAM_THREADS` environment variable
- Async Updates: Non-blocking GUI updates during analysis execution
- Smart Progress Tracking: Streamlined progress indicators with real-time status updates

## Dependencies

### Core Dependencies

- click>=8.0.0: CLI framework
- watchdog>=3.0.0: File system monitoring
- pysam>=0.21.0: BAM file processing
- pandas>=1.3.0: Data manipulation
- numpy>=1.21.0: Numerical computations
- scipy>=1.7.0: Scientific computing
- ruptures>=1.1.0: Change point detection
- tqdm>=4.64.0: Progress bars
- ray[default]>=2.0.0: Distributed computing
- nicegui>=3.0.4: Advanced GUI framework with enhanced performance

### External Dependencies

- bedtools: For region extraction and BED file operations
- samtools: For BAM file manipulation
- R and Rscript: For statistical analysis and classification

### Git Submodules

- nanoDX: NanoDX analysis tools
- hv_rapidCNS2: Rapid CNS analysis tools

## License

This software is provided "as is", and is for research use only.

robin is distributed under a CC BY-NC 4.0 license. See LICENSE for more information. This license does not override any licenses that may be present in the third party tools used by robin.

## Acknowledgments

This tool uses a range of third party tools and applications including:

- [Sturgeon](https://github.com/marcpaga/sturgeon)
- [Rapid-CNS2](https://link.springer.com/article/10.1007/s00401-022-02415-6)
- [Readfish](https://github.com/LooseLab/readfish)
- [cnv_from_bam](https://github.com/adoni5/cnv_from_bam)
- [methylartist](https://github.com/adamewing/methylartist)

We are grateful to the authors of these tools for their work.

We also thank a lot of people who have contributed to these tools including: Graeme Fox, Simon Deacon, Rory Munro, Satrio Wibowo, Thomas Murray, Inswasti Cahyani, Nadine Holmes, Simon Paine, Stuart Smith and many others from outside Nottingham.

We are particularly grateful to Areeba Patel, Felix Sahm and colleagues for their work on Rapid-CNS2.

This list is non-exhaustive and the software is under active development.

In addition we use:

- Click - Python package for creating command line interfaces
- Watchdog - Python library for monitoring file system events
- pysam - Python interface for SAM/BAM files
- Ray - Distributed computing framework
- NiceGUI - Web-based GUI framework

## Next Steps

- Read the [Installation Guide](installation.md) for detailed setup instructions
- Return to the [Home page](../index.md)

# Original Robin

## Basic Usage

ROBIN can be run from the command line:

```bash
robin --threads 2 -r /path/to/reference/hg38.fa -w /path/to/watchfolder /path/to/output
```

This will start the ROBIN interface, which provides a graphical user interface for real-time nanopore data analysis.

A typical watchfolder on a p2i would be:

```bash
/data/<experiment_id>/<sample_id>
```

ROBIN can watch multiple runs at once. So, if you are running two experiments on the same day on a p2i if you give them the same experiment_id, ROBIN will split the results into two runs in the GUI according to the <sample_id>. You do not (in principle...) need to run separate instances of ROBIN.

## Reference Genome Selection

ROBIN requires the use of the human reference genome in hg38/GRCh38 format. For optimal results, we recommend using the `no_alt_analysis_set` version of the reference genome, which excludes alternative loci and decoy sequences. This version is specifically designed for variant calling and analysis.

### Recommended Reference Genome

The recommended reference genome can be downloaded from NCBI:

- [GRCh38 no_alt_analysis_set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)

This version is recommended because:

- It excludes alternative loci (ALT contigs) which can cause false positive variant calls
- It removes decoy sequences that were added to improve alignment
- It maintains compatibility with most analysis tools and databases
- It provides a cleaner reference for variant calling

For more information about human reference genome selection, see [Which human reference genome to use?](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) by Heng Li.

## Bed File Selection

For the CNS tumor analysis pipeline, ROBIN requires a specific bed file to ensure accurate and consistent analysis. The correct bed file to use is:

- `panel_11092024_5kb_pad.bed`

This bed file is specifically designed for CNS tumor analysis and includes:

- Target regions for CNS tumor analysis
- 5kb padding around each target region for improved coverage analysis
- Optimized regions for variant calling and CNV analysis

### Bed File Location

The bed file is included in the ROBIN repository and can be found at:

```
src/robin/resources/panel_11092024_5kb_pad.bed
```

### Using the Bed File

You can specify the bed file in two ways:

1. Using the command line:

```bash
robin --bed_file /path/to/panel_11092024_5kb_pad.bed [other options]
```

2. In your configuration file:

```ini
[options]
bed_file = /path/to/panel_11092024_5kb_pad.bed
```

## Example Workflow

1. **Start a New Analysis**
   ```bash
   robin --threads 2 -r ~/references/hg38.fa -w ~/datasets/new_data /path/to/output
   ```

2. **Monitor Real-time Data**
   - The interface will automatically detect bam files written in the watch folder location.

3. **View Results**
   - Classification results are displayed in real-time
   - Export options are available for further analysis
   - Reports can be generated automatically

## Command Line Options

ROBIN provides a comprehensive set of command line options. You can view all available options using:

```bash
robin --help
```

### Required Options

- `--threads INTEGER`: Number of threads available for processing
- `--reference` or `-r FILE`: Path to the reference genome and index
- `--centreid TEXT`: Provide an identifier to be used in the experiment name field
- `--basecall_config TEXT`: The basecalling configuration to use
- `--experiment_duration INTEGER`: The experiment run time in hours

### Optional Configuration

- `--config` or `-c FILE`: Read option defaults from the specified INI file [default: config.ini]
- `--port INTEGER`: Port for the GUI interface
- `--target_panel` or `-t [rCNS2|AML]`: Select analysis gene panel (default: rCNS2)
- `--watchfolder` or `-w DIRECTORY`: Directory to watch for new data
- `--bed_file` or `-b FILE`: Path to the bedfile to use for adaptive sampling
- `--kit TEXT`: Specify sequencing kit
- `--force_sampleid TEXT`: Force a specific sampleID
- `--readfish_toml` or `-rt FILE`: Path to the TOML file used to control readfish

### Analysis Control

- `--exclude` or `-e [sturgeon|forest|nanodx|pannanodx|cnv|fusion|coverage|mgmt]`: Exclude specific analysis types
- `--enable-snp-calling`: Enable SNP calling functionality (requires reference genome)
- `--sequencing_summary FILE`: Path to sequencing summary file for timestamp information
- `--simtime BOOLEAN`: Simulate the addition of existing files to the pipeline based on read data

### Logging and Debugging

- `--log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL]`: Set the logging level
- `--log-file FILE`: Path to the log file [default: ROBIN.log]
- `--showerrors / --noerrors`: Control display of R errors
- `--browse`: Browse Historic Data
- `--no-telemetry`: Opt out of sending anonymous usage statistics

### Other Options

- `--version`: Show the version and exit
- `--help`: Show help message and exit

## Configuration

ROBIN can be configured using an INI file. A basic configuration might look like:

```ini
[options]
# ROBIN needs to be aware of the targets being used in your experiment
bed_file = /path/to/bed/file/
# The centreID is displayed in PDF reports generated by ROBIN
centreID = NUH
# This is the port that the GUI will be displayed
port = 5678
# ROBIN needs to know where the reference file for aligning is
reference = /path/to/reference/hg38.fasta
# You can vary the number of threads used by ROBIN. For P2i devices this should be 1
threads = 8

# The following options will be used in future improvements to ROBIN

# We record the basecall configuration in the options for future features
basecall_config = dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg
# We track the ideal experiment duration for future features
experiment_duration = 24
# We log the experiment kit being used for future features
kit = SQK-RAD114
```

To use a config file, run ROBIN with the config argument:

```bash
robin -c config.ini -w /path/to/watchfolder /path/to/output
```

### Output Directory

The output directory (`/path/to/output`) is where ROBIN will store all analysis results. Before running ROBIN:

1. Create this directory if it doesn't exist
2. Ensure you have write permissions to this location

ROBIN will automatically create subfolders within this directory, with each subfolder named according to the sampleID assigned in MinKNOW. For example, if your MinKNOW sampleID is "CNS_2024_001", ROBIN will create a folder structure like:

```
/path/to/output/
|-- CNS_2024_001/
|   `-- real-time analysis results
```

For example:

```bash
# Create the main output directory
mkdir -p ~/robin_analysis
# ROBIN will create sampleID-named subfolders automatically based on MinKNOW
robin -c config.ini -w /path/to/watchfolder ~/robin_analysis
```

## Important Notes

- The reference genome must be in hg38 format
- For P2i devices, set threads to 1
- The watchfolder is where minKNOW will be writing the aligned BAM files
- The output folder must exist and should be empty

## Next Steps

- Read the [Installation Guide](installation.md) for detailed setup instructions
- Check the [GitHub repository](https://github.com/LooseLab/robin) for updates and issues
- Contact us via [this form](https://forms.gle/kdX2eiPQPdDUpaBE9) for support or collaboration
# ROBIN Quickstart

This guide covers the essential steps to run ROBIN (Little John version).

## Usage

ROBIN expects to analyse BAM files generated during sequencing by an Oxford Nanopore Technologies sequencer. ROBIN presumes real time HAC basecalling (SUP is not required). ROBIN presumes data have been called with 5hmC 5mC methylation calling in MinKNOW. ROBIN presumes real time alignment is running in MinKNOW - ROBIN does not realign your reads.

ROBIN assumes that BAM files are being output in small batches. We recommend setting file output to one bam for every 10,000 to 50,000 reads. We do not support real time processing of BAM files in 1 hour chunks (the default output).

ROBIN does not use pod5 data or fastq data from the sequencer - you can deselect these options if you wish.

## Important

On platforms with 64Gb of RAM or less we recommend restarting your device prior to a run. As an example, if you are running a p2i and start a run on position A and later on position B we recommend you restart at the end of the run on position B (once base calling is complete). You can simply restart dorado if you know how to do this. We find dorado holds on to memory for an indefinite period of time and this can cause problems.

## Run a Workflow

The primary command for running robin workflows is:

```bash
robin workflow <data_folder> --work-dir <output_folder> -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest --reference ~/references/hg38_simple.fa --center <center_id>
```

### Command Breakdown

- `robin workflow`: The main workflow command
- `<data_folder>`: Directory containing your BAM files
- `--work-dir <output_folder>`: Directory where results will be saved
- `-w`: Workflow specification (comma-separated list of analysis types)
- `--reference`: Path to reference genome (required for some analyses)
- `--center <center_id>`: Center ID running the analysis (e.g., `Sherwood`, `Auckland`, `New York`)
- `--target-panel`: The specific panel that is being applied

## Example Usage

```bash
# Basic workflow with all analysis types
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion,mgmt,sturgeon,nanodx,pannanodx,random_forest \
  --reference ~/references/hg38_simple.fa \
  --center Sherwood \
  --target-panel rCNS2

# Simplified workflow with just a few analyses
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center Auckland \
  --target-panel PanCan

# With verbose output and custom logging
robin workflow ~/data/bam_files \
  --work-dir ~/results \
  -w mgmt,cnv,sturgeon \
  --reference ~/references/hg38_simple.fa \
  --center New_York \
  --target-panel rCNS2 \
  --verbose \
  --log-level INFO
```

## Known Issues

- Currently SNP calling is not enabled in this version of ROBIN. It will be re-enabled in the near future.
- CNV change inference is based on extensive heuristics - every call should be checked by visual inspection.
- If you ctrl-c to end ROBIN it will do its best to clean up and stop gracefully but may fail.
- CSV data export is in development but is not currently available - it will be enabled in the near future.
- Many other unknown issues - please open an issue and we will resolve where possible.

## Available Commands

### list-job-types

List all available job types organized by queue category:

```bash
robin list-job-types
```

Available Job Types:

- Preprocessing Queue: preprocessing
- Bed Conversion Queue: bed_conversion
- Analysis Queue: mgmt, cnv, target, fusion
- Classification Queue: sturgeon, nanodx, pannanodx
- Slow Queue: random_forest

### workflow

Run an async workflow on BAM files in a directory:

```bash
robin workflow /path/to/directory --workflow "workflow_plan" [OPTIONS]
```

Required Options:

- `--workflow`, `-w`: Workflow plan (e.g., `mgmt,sturgeon` or `preprocessing:bed_conversion,analysis:mgmt,classification:sturgeon`)
- `--center`: Center ID running the analysis (e.g., `Sherwood`, `Auckland`, `New York`)

Optional Options:

- `--work-dir`, `-d`: Base output directory for analysis results
- `--reference`, `-r`: Path to reference genome (FASTA format)
- `--verbose`, `-v`: Enable verbose output and detailed error traces
- `--no-process-existing`: Skip processing existing files, only watch for new changes
- `--log-level`: Global log level (DEBUG|INFO|WARNING|ERROR, default: ERROR)
- `--job-log-level`: Set log level for specific job (e.g., `preprocessing:DEBUG`, `mgmt:WARNING`)
- `--deduplicate-jobs`: Job types to deduplicate by sample ID (e.g., `sturgeon`, `mgmt`)
- `--no-progress`: Disable progress bars for file processing
- `--use-ray`/`--no-use-ray`: Enable Ray distributed computing (default: on)
- `--with-gui`/`--no-gui`: Launch NiceGUI workflow monitor (default: on)

## Panel Management

Manage built-in and custom target panels used by analyses like target, cnv, and fusion.

Built-in panels include: rCNS2, AML, PanCan.
Custom panels are stored internally after you add them from a BED file.

### List available panels

```bash
robin list-panels
```

### Add a custom panel from a BED file

```bash
# Add and register a panel (BED must have >=4 columns: chr, start, end, gene_name[s])
robin add-panel /path/to/your_panel.bed MyCustomPanel

# Optional: validate format only, without adding
robin add-panel /path/to/your_panel.bed MyCustomPanel --validate-only
```

Notes:

- Panel names cannot be empty and cannot reuse reserved names: rCNS2, AML, PanCan.
- BED may be 4- or 6-column; if multiple genes are in one region, use comma-separated names.

### Remove a custom panel

```bash
# Will prompt for confirmation
robin remove-panel MyCustomPanel

# Skip confirmation
robin remove-panel MyCustomPanel --force
```

Built-in panels cannot be removed.

### Use a panel in a workflow

```bash
robin workflow /path/to/bam_files \
  --work-dir ~/results \
  -w target,cnv,fusion \
  --target-panel MyCustomPanel \
  --center Sherwood
```

## Performance Features

### Enhanced Processing Capabilities

- Batched Processing: All analysis workflows now support batched processing for improved efficiency
- Memory Optimization: Intelligent memory management for large datasets and long-running processes
- Multi-threading: Configurable multi-threaded BAM processing via `LJ_BAM_THREADS` environment variable
- Async Updates: Non-blocking GUI updates during analysis execution
- Smart Progress Tracking: Streamlined progress indicators with real-time status updates

## Dependencies

### Core Dependencies

- click>=8.0.0: CLI framework
- watchdog>=3.0.0: File system monitoring
- pysam>=0.21.0: BAM file processing
- pandas>=1.3.0: Data manipulation
- numpy>=1.21.0: Numerical computations
- scipy>=1.7.0: Scientific computing
- ruptures>=1.1.0: Change point detection
- tqdm>=4.64.0: Progress bars
- ray[default]>=2.0.0: Distributed computing
- nicegui>=3.0.4: Advanced GUI framework with enhanced performance

### External Dependencies

- bedtools: For region extraction and BED file operations
- samtools: For BAM file manipulation
- R and Rscript: For statistical analysis and classification

### Git Submodules

- nanoDX: NanoDX analysis tools
- hv_rapidCNS2: Rapid CNS analysis tools

## License

This software is provided "as is", and is for research use only.

robin is distributed under a CC BY-NC 4.0 license. See LICENSE for more information. This license does not override any licenses that may be present in the third party tools used by robin.

## Acknowledgments

This tool uses a range of third party tools and applications including:

- [Sturgeon](https://github.com/marcpaga/sturgeon)
- [Rapid-CNS2](https://link.springer.com/article/10.1007/s00401-022-02415-6)
- [Readfish](https://github.com/LooseLab/readfish)
- [cnv_from_bam](https://github.com/adoni5/cnv_from_bam)
- [methylartist](https://github.com/adamewing/methylartist)

We are grateful to the authors of these tools for their work.

We also thank a lot of people who have contributed to these tools including: Graeme Fox, Simon Deacon, Rory Munro, Satrio Wibowo, Thomas Murray, Inswasti Cahyani, Nadine Holmes, Simon Paine, Stuart Smith and many others from outside Nottingham.

We are particularly grateful to Areeba Patel, Felix Sahm and colleagues for their work on Rapid-CNS2.

This list is non-exhaustive and the software is under active development.

In addition we use:

- Click - Python package for creating command line interfaces
- Watchdog - Python library for monitoring file system events
- pysam - Python interface for SAM/BAM files
- Ray - Distributed computing framework
- NiceGUI - Web-based GUI framework

## Next Steps

- Read the [Installation Guide](installation.md) for detailed setup instructions
- Return to the [Home page](../index.md)
# ROBIN Quickstart

This guide will help you get started with ROBIN quickly.

## Basic Usage

ROBIN can be run from the command line:

```bash
robin --threads 2 -r /path/to/reference/hg38.fa -w /path/to/watchfolder /path/to/output
```

This will start the ROBIN interface, which provides a graphical user interface for real-time nanopore data analysis.

A typical watchfolder on a p2i would be:

```bash
/data/<experiment_id>/<sample_id>
```

ROBIN can watch multiple runs at once. So, if you are running two experiments on the same day on a p2i if you give them the same experiment_id, ROBIN will split the results into two runs in the GUI according to the <sample_id>. You do not (in principle...) need to run separate instances of ROBIN.

## Reference Genome Selection

ROBIN requires the use of the human reference genome in hg38/GRCh38 format. For optimal results, we recommend using the `no_alt_analysis_set` version of the reference genome, which excludes alternative loci and decoy sequences. This version is specifically designed for variant calling and analysis.

### Recommended Reference Genome

The recommended reference genome can be downloaded from NCBI:
- [GRCh38 no_alt_analysis_set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)

This version is recommended because:
- It excludes alternative loci (ALT contigs) which can cause false positive variant calls
- It removes decoy sequences that were added to improve alignment
- It maintains compatibility with most analysis tools and databases
- It provides a cleaner reference for variant calling

For more information about human reference genome selection, see [Which human reference genome to use?](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) by Heng Li.

## Bed File Selection

For the CNS tumor analysis pipeline, ROBIN requires a specific bed file to ensure accurate and consistent analysis. The correct bed file to use is:

- `panel_11092024_5kb_pad.bed`

This bed file is specifically designed for CNS tumor analysis and includes:
- Target regions for CNS tumor analysis
- 5kb padding around each target region for improved coverage analysis
- Optimized regions for variant calling and CNV analysis

### Bed File Location

The bed file is included in the ROBIN repository and can be found at:
```
src/robin/resources/panel_11092024_5kb_pad.bed
```


### Using the Bed File

You can specify the bed file in two ways:

1. Using the command line:
```bash
robin --bed_file /path/to/panel_11092024_5kb_pad.bed [other options]
```

2. In your configuration file:
```ini
[options]
bed_file = /path/to/panel_11092024_5kb_pad.bed
```

## Example Workflow

1. **Start a New Analysis**
   ```bash
   robin --threads 2 -r ~/references/hg38.fa -w ~/datasets/new_data /path/to/output
   ```

2. **Monitor Real-time Data**
   - The interface will automatically detect bam files written in the watch folder location.
   
3. **View Results**
   - Classification results are displayed in real-time
   - Export options are available for further analysis
   - Reports can be generated automatically

## Command Line Options

ROBIN provides a comprehensive set of command line options. You can view all available options using:

```bash
robin --help
```

### Required Options

- `--threads INTEGER`: Number of threads available for processing
- `--reference` or `-r FILE`: Path to the reference genome and index
- `--centreid TEXT`: Provide an identifier to be used in the experiment name field
- `--basecall_config TEXT`: The basecalling configuration to use
- `--experiment_duration INTEGER`: The experiment run time in hours

### Optional Configuration

- `--config` or `-c FILE`: Read option defaults from the specified INI file [default: config.ini]
- `--port INTEGER`: Port for the GUI interface
- `--target_panel` or `-t [rCNS2|AML]`: Select analysis gene panel (default: rCNS2)
- `--watchfolder` or `-w DIRECTORY`: Directory to watch for new data
- `--bed_file` or `-b FILE`: Path to the bedfile to use for adaptive sampling
- `--kit TEXT`: Specify sequencing kit
- `--force_sampleid TEXT`: Force a specific sampleID
- `--readfish_toml` or `-rt FILE`: Path to the TOML file used to control readfish

### Analysis Control

- `--exclude` or `-e [sturgeon|forest|nanodx|pannanodx|cnv|fusion|coverage|mgmt]`: Exclude specific analysis types
- `--enable-snp-calling`: Enable SNP calling functionality (requires reference genome)
- `--sequencing_summary FILE`: Path to sequencing summary file for timestamp information
- `--simtime BOOLEAN`: Simulate the addition of existing files to the pipeline based on read data

### Logging and Debugging

- `--log-level [DEBUG|INFO|WARNING|ERROR|CRITICAL]`: Set the logging level
- `--log-file FILE`: Path to the log file [default: ROBIN.log]
- `--showerrors / --noerrors`: Control display of R errors
- `--browse`: Browse Historic Data
- `--no-telemetry`: Opt out of sending anonymous usage statistics

### Other Options

- `--version`: Show the version and exit
- `--help`: Show help message and exit



## Configuration

ROBIN can be configured using an INI file. A basic configuration might look like:

```ini
[options]
# ROBIN needs to be aware of the targets being used in your experiment
bed_file = /path/to/bed/file/
# The centreID is displayed in PDF reports generated by ROBIN
centreID = NUH
# This is the port that the GUI will be displayed
port = 5678
# ROBIN needs to know where the reference file for aligning is
reference = /path/to/reference/hg38.fasta
# You can vary the number of threads used by ROBIN. For P2i devices this should be 1
threads = 8

# The following options will be used in future improvements to ROBIN

# We record the basecall configuration in the options for future features
basecall_config = dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg
# We track the ideal experiment duration for future features
experiment_duration = 24
# We log the experiment kit being used for future features
kit = SQK-RAD114
```

To use a config file, run ROBIN with the config argument:

```bash
robin -c config.ini -w /path/to/watchfolder /path/to/output
```

### Output Directory

The output directory (`/path/to/output`) is where ROBIN will store all analysis results. Before running ROBIN:

1. Create this directory if it doesn't exist
2. Ensure you have write permissions to this location

ROBIN will automatically create subfolders within this directory, with each subfolder named according to the sampleID assigned in MinKNOW. For example, if your MinKNOW sampleID is "CNS_2024_001", ROBIN will create a folder structure like:

```
/path/to/output/
??? CNS_2024_001/
    ??? real-time analysis results
```

For example:
```bash
# Create the main output directory
mkdir -p ~/robin_analysis
# ROBIN will create sampleID-named subfolders automatically based on MinKNOW
robin -c config.ini -w /path/to/watchfolder ~/robin_analysis
```

## Important Notes

- The reference genome must be in hg38 format
- For P2i devices, set threads to 1
- The watchfolder is where minKNOW will be writing the aligned BAM files
- The output folder must exist and should be empty

## Next Steps

- Read the [Installation Guide](installation.md) for detailed setup instructions
- Check the [GitHub repository](https://github.com/LooseLab/robin) for updates and issues
- Contact us via [this form](https://forms.gle/kdX2eiPQPdDUpaBE9) for support or collaboration 