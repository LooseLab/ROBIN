# Command-line reference

This section documents the `robin` command-line interface. For a first run, use the [Quickstart](../getting-started/quickstart.md); for installation, see [Installation](../getting-started/installation.md).

Run `robin --help` and `robin <command> --help` to confirm options in your installed version.

## Commands

<div class="grid" markdown>

<div class="robin-feature-card" markdown>
### [`robin workflow`](workflow.md)
Watch a BAM directory and run the selected pipeline using Ray or threaded execution, with optional NiceGUI monitoring.
</div>

<div class="robin-feature-card" markdown>
### [`robin list-job-types`](jobs.md)
Show available analysis job types and how they map to workflow queues.
</div>

<div class="robin-feature-card" markdown>
### [Panel commands](panels.md)
List built-in panels and add or remove custom BED panels.
</div>

<div class="robin-feature-card" markdown>
### [`robin utils`](utils.md)
Stage references and panels, update model/ClinVar resources and run utility analyses.
</div>

<div class="robin-feature-card" markdown>
### [`robin password`](password.md)
Set or replace the default administrator password used by the web interface.
</div>

</div>

## Startup and consent

Commands that start protected processing may display the research-use disclaimer and require `I agree`. The workflow can also prompt for initial GUI-user setup.

See [Starting ROBIN](../getting-started/startup.md) for the complete startup sequence.

## Selected environment variables

| Variable | Effect |
|----------|--------|
| `ROBIN_PROCESS_LARGE_BAMS` | Enables large-BAM behaviour intended for non-live processing; ROBIN warns against combining this with live sequencing. |
| `LJ_BAM_THREADS` | Controls optional BAM decompression/read threading. |

Environment variables specific to individual analyses are documented with those analyses or their job configuration.

## Web interface

For browser navigation, authentication, sample pages and result interpretation, see [Using ROBIN](../using-robin/index.md).

## Related

- [`robin workflow`](workflow.md)
- [Job types](jobs.md)
- [Panel commands](panels.md)
- [`robin utils`](utils.md)
