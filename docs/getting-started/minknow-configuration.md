# MinKNOW Configuration Guide

This guide explains how to configure MinKNOW for use with ROBIN, including settings for adaptive sampling and real-time analysis.

## Overview

ROBIN requires specific MinKNOW settings to enable:

- Real-time basecalling

- Adaptive sampling

- BAM file output

- Proper sample identification

## Required MinKNOW Settings

### 1. Basecalling Configuration

ROBIN has been designed to work with high accuracy basecalling and requires methyaltion data. MinKNOW provides 5hmc and 5mc cg calling.

*IMPORTANT* ROBIN uses methylation calling data from off target reads. You must not switch off basecalling for these reads.

At the time of writing, the following model is the one used for base calling.
- `dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_prom.cfg`


### 2. Adaptive Sampling Settings

Configure adaptive sampling with the following settings:

- Use the bed file: `panel_11092024_5kb_pad.bed`

- Use the same reference file as used by ROBIN.

- Use the "enrich" mode.


### 3. Alignment Settings

ROBIN requires BAM files to be aligned to a reference.

- Use the same reference file as for the Adaptive Sampling Settings and ROBIN.


### 4. Output Settings

Configure the following output settings:

- Ensure BAM file output and set to either 8000 reads per file or one file per minute.

- Optional - disable pod5 output. 

- Disable FastQ output.

### 5. Sample Settings

For each run:

- Set a unique sample ID (this will be used by ROBIN for output organization)

- Any experiment ID can be used. ROBIN ignores this. However, it can be useful to use the experiment ID to group related experiments and use this folder to point ROBIN too to monitor specific subsets of data.

= Set the experiment duration (typically 24 hours for CNS tumor analysis)


## Starting a Run

- Configure MinKNOW with the settings above

- Begin the MinKNOW run

- Start ROBIN and give it the output folder that MinKNOW will be writing data to.

- ROBIN will automatically detect and process the BAM files

## Troubleshooting

## Next Steps

After configuring MinKNOW:
1. Start ROBIN as described in the [ROBIN Quickstart](quickstart.md)

1. Begin your MinKNOW run

1. Monitor the analysis in the ROBIN interface

## Additional Resources

- [MinKNOW Documentation](https://nanoporetech.com/minknow)
- [Adaptive Sampling Guide](https://nanoporetech.com/adaptive-sampling)
- [ROBIN Quickstart](quickstart.md) 