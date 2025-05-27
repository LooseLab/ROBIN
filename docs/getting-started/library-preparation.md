# Library Preparation Guide

This guide provides detailed instructions for preparing libraries for analysis with ROBIN.

## Overview

ROBIN is designed to work with libraries prepared using the Oxford Nanopore Technologies (ONT) sequencing platform. 

For CNS tumor analysis, you can use either:

1. The ROBIN-specific library preparation protocol
1. The standard ONT LSK114 ligation protocol

## Input DNA

As ROBIN uses adaptive sampling, we aim for read lengths in the range of 8-12 kb.

## Library Preparation Options

### Option 1: ROBIN Protocol

The ROBIN-specific library preparation protocol is detailed in our protocols.io publication:
[ROBIN: A Unified Nanopore-based Sequencing Assay](https://www.protocols.io/view/robin-a-unified-nanopore-based-sequencing-assay-in-bp2l6xepklqe/v3)

This protocol has been optimized for CNS tumor analysis and includes:

- Specific fragmentation conditions
- Optimized adapter ligation steps
- Quality control checkpoints
- Detailed troubleshooting guidance

### Option 2: Standard ONT Protocol

Alternatively, you can use the standard ONT LSK114 ligation protocol, which is available from:
- [ONT LSK114 Protocol Documentation](https://store.nanoporetech.com/uk/ligation-sequencing-kit-v14.html)

## DNA Input Requirements

- Input DNA should be of high quality (A260/A280 ratio between 1.8-2.0)
- Minimum input: 100ng
- Optimal input: 500ng-1μg
- DNA should be in TE buffer or nuclease-free water



## Next Steps

After library preparation:

1. Proceed to sequencing setup and configure MinKNOW for your run

1. Start ROBIN analysis as described in the [ROBIN Quickstart](quickstart.md)

## Additional Resources

- [ROBIN Protocol on protocols.io](https://www.protocols.io/view/robin-a-unified-nanopore-based-sequencing-assay-in-bp2l6xepklqe/v3)
- [ONT Library Preparation Documentation](https://nanoporetech.com/library-preparation)
- [SQK-RAD114 Protocol](https://store.nanoporetech.com/sqk-rad114.html)
- [ROBIN Quickstart](quickstart.md) 