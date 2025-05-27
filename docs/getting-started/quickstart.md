# ROBIN Quickstart

This guide will help you get started with ROBIN quickly.

## Basic Usage

ROBIN can be run from the command line:

```bash
robin --threads 2 -r /path/to/reference/hg38.fa -w /path/to/watchfolder /path/to/output
```

This will start the ROBIN interface, which provides a graphical user interface for real-time nanopore data analysis.

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
└── CNS_2024_001/
    └── real-time analysis results
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