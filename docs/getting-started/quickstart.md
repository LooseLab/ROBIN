# Quick Start Guide

This guide will help you get started with ROBIN quickly.

## Basic Usage

ROBIN can be run from the command line:

```bash
robin run
```

This will start the ROBIN interface, which provides a graphical user interface for real-time nanopore data analysis.

## Command Line Options

ROBIN provides several command line options:

```bash
robin --help
```

Common options include:

- `--config`: Specify a configuration file
- `--output`: Set the output directory
- `--device`: Specify the sequencing device
- `--model`: Select the classification model

## Example Workflow

1. **Start a New Analysis**
   ```bash
   robin run --config my_config.yaml
   ```

2. **Monitor Real-time Data**
   - The interface will automatically connect to your sequencing device
   - Real-time classification results will be displayed
   - Data is automatically saved for later analysis

3. **View Results**
   - Classification results are displayed in real-time
   - Export options are available for further analysis
   - Reports can be generated automatically

## Configuration

ROBIN uses YAML configuration files. A basic configuration might look like:

```yaml
device:
  type: "minion"
  host: "localhost"
  port: 8000

analysis:
  model: "default"
  batch_size: 1000
  confidence_threshold: 0.8

output:
  directory: "results"
  format: "csv"
```

For more detailed configuration options, see the [Configuration Guide](../user-guide/configuration.md).

## Next Steps

- Read the [User Guide](../user-guide/overview.md) for detailed information
- Explore the [API Reference](../api/core.md) for programmatic usage
- Check the [Configuration Guide](../user-guide/configuration.md) for advanced settings 